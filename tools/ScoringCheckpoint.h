#ifndef BRIQ_SCORING_CHECKPOINT_H
#define BRIQ_SCORING_CHECKPOINT_H

#include <cstdint>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <fcntl.h>
#include <unistd.h>

// Non-cryptographic integrity check for accidental input changes on resume.
// Hash paths and complete file contents; never infer identity from filenames.
static std::string batchFingerprint(const std::vector<std::string>& jobs) {
    uint64_t hash=14695981039346656037ULL;
    auto add=[&](unsigned char byte){hash^=byte;hash*=1099511628211ULL;};
    for (size_t i=0;i<jobs.size();++i) {
        for (size_t j=0;j<jobs[i].size();++j) add(jobs[i][j]);
        add(0);
        std::ifstream file(jobs[i].c_str(),std::ios::binary);
        if (!file) throw std::runtime_error("Cannot fingerprint batch input");
        char buffer[16384];
        while (file.read(buffer,sizeof(buffer)) || file.gcount())
            for (std::streamsize j=0;j<file.gcount();++j) add(buffer[j]);
        if (!file.eof()) throw std::runtime_error("Failed reading batch input fingerprint");
        add(0);
    }
    std::ostringstream out;out<<std::hex<<std::setw(16)<<std::setfill('0')<<hash;
    return out.str();
}

static std::vector<std::string> readCsvRow(const std::string& line) {
    std::vector<std::string> result;
    size_t i=0;
    while (true) {
        std::string value;
        if (i<line.size() && line[i]=='"') {
            ++i;bool closed=false;
            while (i<line.size()) {
                char c=line[i++];
                if (c=='"') {
                    if (i<line.size() && line[i]=='"') {value+='"';++i;}
                    else {closed=true;break;}
                } else value+=c;
            }
            if (!closed) throw std::runtime_error("Truncated quoted CSV field");
        } else {
            while (i<line.size() && line[i]!=',') {
                if (line[i]=='"') throw std::runtime_error("Invalid CSV quoting");
                value+=line[i++];
            }
        }
        result.push_back(value);
        if (i==line.size()) break;
        if (line[i++]!=',') throw std::runtime_error("Invalid CSV separator");
    }
    return result;
}

static std::vector<double> readCompletedTotals(const std::string& path,
                                              const std::vector<std::string>& jobs) {
    std::ifstream file(path.c_str(),std::ios::binary);
    std::string line;
    if (!std::getline(file,line)) throw std::runtime_error("Missing totals CSV header");
    if (!line.empty() && line.back()=='\r') line.pop_back();
    if (line!="job_index,input_pdb,E_total") throw std::runtime_error("Invalid totals CSV header");
    std::vector<double> energies;
    while (std::getline(file,line)) {
        if (!line.empty() && line.back()=='\r') line.pop_back();
        std::vector<std::string> fields=readCsvRow(line);
        size_t index=energies.size();
        std::ostringstream expected;expected<<index;
        double energy;
        if (fields.size()!=3 || index>=jobs.size() || fields[0]!=expected.str() ||
            fields[1]!=jobs[index] || !parseDouble(fields[2].c_str(),energy))
            throw std::runtime_error("Totals CSV is not a finite, ordered prefix of this manifest");
        energies.push_back(energy);
    }
    if (!file.eof()) throw std::runtime_error("Cannot read existing totals CSV");
    file.clear();file.seekg(-1,std::ios::end);char last=0;file.get(last);
    if (last!='\n') throw std::runtime_error("Totals CSV has an incomplete final line; preserve and inspect it before resuming");
    return energies;
}

static std::string checkpointText(const std::vector<std::string>& jobs,
                                  const std::string& engine) {
    std::ostringstream out;
    out<<"BRiQ-total-checkpoint-v1\nsource="<<BRIQ_SOURCE_VERSION
       <<"\nengine="<<engine<<"\nseed=0\nweights=1,1,1\ncount="<<jobs.size()
       <<"\ninputs-fnv1a64="<<batchFingerprint(jobs)<<"\n";
    return out.str();
}

static void createCheckpoint(const std::string& path,const std::string& text) {
    int descriptor=open(path.c_str(),O_WRONLY|O_CREAT|O_EXCL,0644);
    if (descriptor<0) throw std::runtime_error("Cannot exclusively create checkpoint: "+path);
    size_t offset=0;
    while (offset<text.size()) {
        ssize_t written=write(descriptor,text.data()+offset,text.size()-offset);
        if (written<=0) {close(descriptor);throw std::runtime_error("Cannot write checkpoint");}
        offset+=written;
    }
    if (close(descriptor)!=0) throw std::runtime_error("Cannot close checkpoint");
}

static uint64_t csvHashUpdate(uint64_t hash,const std::string& text) {
    for (size_t i=0;i<text.size();++i) {hash^=static_cast<unsigned char>(text[i]);hash*=1099511628211ULL;}
    return hash;
}
static std::vector<uint64_t> csvPrefixHashes(const std::string& path) {
    std::ifstream file(path.c_str(),std::ios::binary);
    std::string line;
    std::vector<uint64_t> hashes;
    uint64_t hash=14695981039346656037ULL;
    while (std::getline(file,line)) {hash=csvHashUpdate(hash,line+"\n");hashes.push_back(hash);}
    if (hashes.empty()) throw std::runtime_error("Cannot hash totals CSV");
    return hashes;
}
static void recordCheckpoint(const std::string& path,size_t count,uint64_t hash) {
    std::ofstream file(path.c_str(),std::ios::app);
    file<<count<<" "<<std::hex<<hash<<"\n";
    file.close();
    if (!file.good()) throw std::runtime_error("Cannot record committed CSV prefix");
}
static size_t verifiedCheckpointPrefix(const std::string& content,
        const std::string& identity,const std::vector<uint64_t>& hashes) {
    if (content.compare(0,identity.size(),identity)!=0)
        throw std::runtime_error("Checkpoint source, engine or input contents do not match");
    std::istringstream records(content.substr(identity.size()));
    std::string line;
    size_t committed=0;
    bool any=false;
    while (std::getline(records,line)) {
        std::istringstream fields(line);
        size_t count;uint64_t hash;std::string extra;
        if (!(fields>>count>>std::hex>>hash) || (fields>>extra) ||
            count>=hashes.size() || (any && count<=committed) || hashes[count]!=hash)
            throw std::runtime_error("Invalid checkpoint or modified totals CSV prefix");
        committed=count;any=true;
    }
    if (!content.empty() && content.back()!='\n')
        throw std::runtime_error("Incomplete checkpoint record; preserve and inspect it before resuming");
    return committed;
}
#endif
