#ifndef BRIQ_MUTATION_BATCH_H
#define BRIQ_MUTATION_BATCH_H

#include "model/RNABaseLib.h"
#include <map>
#include <tuple>

// Included after the common CLI/scoring helpers in RunEnergy.cpp.
namespace {
using MutationKey = tuple<string, string, char>;
using ResidueIdentity = pair<string, string>;

static vector<string> mutationFields(string line) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    vector<string> fields;
    size_t start = 0;
    while (true) {
        size_t end = line.find('\t', start);
        fields.push_back(line.substr(start, end == string::npos ? end : end-start));
        if (end == string::npos) return fields;
        start = end+1;
    }
}

static string absoluteMutationPath(const string& path) {
    char resolved[4096];
    if (realpath(path.c_str(), resolved) == NULL)
        throw runtime_error("Cannot resolve mutation input: " + path);
    return resolved;
}

static string mutationChain(const string& name) { return name.empty() ? " " : name; }

struct MutationGeometry {
    map<MutationKey, unique_ptr<RNABase> > bases;
    explicit MutationGeometry(const string& path) {
        ifstream input(path.c_str());
        string line;
        if (!getline(input, line) || mutationFields(line) != vector<string>(
                {"chain", "residue", "base", "atom", "x", "y", "z"}))
            throw runtime_error("Invalid geometry TSV header");
        size_t number = 1;
        while (getline(input, line)) {
            ++number;
            vector<string> f = mutationFields(line);
            if (f.size()!=7 || f[0].size()>1 || f[1].empty() || f[2].size()!=1 ||
                    string("ACGU").find(f[2][0])==string::npos || f[3].empty())
                throw runtime_error("Invalid geometry TSV row " + to_string(number));
            double x,y,z;
            if (!parseDouble(f[4].c_str(),x) || !parseDouble(f[5].c_str(),y) ||
                    !parseDouble(f[6].c_str(),z))
                throw runtime_error("Non-finite/invalid geometry coordinate");
            MutationKey key(mutationChain(f[0]),f[1],f[2][0]);
            auto& base = bases[key];
            if (!base) base.reset(new RNABase(f[1],mutationChain(f[0]),f[2][0]));
            if (base->hasAtom(f[3])) throw runtime_error("Duplicate cached atom");
            // Match Atom(PDB-line), which parses coordinates through float.
            // Extra precision here can cross statistical-table bin boundaries.
            if (fabs(x)>numeric_limits<float>::max() || fabs(y)>numeric_limits<float>::max() ||
                    fabs(z)>numeric_limits<float>::max()) throw runtime_error("Coordinate outside PDB numeric range");
            base->addAtom(new Atom(f[3],XYZ(static_cast<float>(x),static_cast<float>(y),static_cast<float>(z))));
        }
        if (!input.eof() || bases.empty()) throw runtime_error("Empty/unreadable geometry cache");
        RnaAtomLib atoms;
        for (auto& entry : bases) {
            RNABase* base = entry.second.get();
            if (!base->hasAtom("C1'") || !base->sidechainComplete(&atoms))
                throw runtime_error("Cached residue lacks base/frame atoms");
            base->updateCoordSystem();
            if (!base->hasLocalFrame) throw runtime_error("Cached residue has no base frame");
            // The file scorer standardizes raw PDB atoms before computing the
            // frame. Validation must not leave the pre-standardization frame
            // cached: RNABaseLib deliberately preserves hasLocalFrame.
            base->hasLocalFrame=false;
        }
    }
};

struct CachedMutationNode {
    // Node does not own its rotamer. Destruction order matters.
    unique_ptr<BaseRotamer> rotamer;
    unique_ptr<BRNode> node;
};

// Always restore original pointers before BRFoldingTree destroys its nodes.
class MutationNodeGuard {
    BRFoldingTree& tree;
    vector<BRNode*> original;
 public:
    explicit MutationNodeGuard(BRFoldingTree& t): tree(t), original(t.nodes,t.nodes+t.seqLen) {}
    void restore() {
        for (int i=0;i<tree.seqLen;++i) {
            tree.nodes[i]=original[i];
            tree.seq[i]=original[i]->baseType;
        }
    }
    ~MutationNodeGuard() { restore(); }
};

class MutationWorkingDirectory {
    string previous;
 public:
    MutationWorkingDirectory() {
        char buffer[4096];
        if (!getcwd(buffer,sizeof(buffer))) throw runtime_error("Cannot read working directory");
        previous=buffer;
    }
    ~MutationWorkingDirectory() { if (chdir(previous.c_str())!=0) cerr<<"Cannot restore cwd"<<endl; }
};

static void validateMutationTopology(const vector<RNABase*>& reference,
                                     const MutationGeometry& geometry) {
    map<ResidueIdentity,size_t> identities;
    vector<vector<RNABase*> > choices(reference.size());
    for (size_t i=0;i<reference.size();++i) {
        ResidueIdentity id(reference[i]->chainID,reference[i]->baseID);
        if (!identities.emplace(id,i).second) throw runtime_error("Ambiguous reference residue identifier");
        choices[i].push_back(reference[i]);
    }
    for (const auto& entry : geometry.bases) {
        auto found=identities.find(ResidueIdentity(get<0>(entry.first),get<1>(entry.first)));
        if (found==identities.end()) throw runtime_error("Cached residue not found in reference");
        choices[found->second].push_back(entry.second.get());
    }
    // The fast path retains reference sequence separation. Refuse any geometry
    // that would change breaks under the existing PDB-scoring preparation.
    for (size_t i=0;i+1<reference.size();++i) {
        if (reference[i]->chainID!=reference[i+1]->chainID) continue;
        const bool connected=reference[i]->connectToNeighbor(*reference[i+1]);
        for (RNABase* a : choices[i]) for (RNABase* b : choices[i+1])
            if (a->connectToNeighbor(*b)!=connected)
                throw runtime_error("Cached geometry changes reference chain connectivity");
    }
}

// Exact memoization: phosphate building depends only on these immutable node
// geometries (cs2, ribose coordinates, both improper angles), not on other bases.
// Cache the result for each ordered pair of reference/cached node identities.
class MutationPhosphates {
    map<pair<int,int>, PhophateGroupLocal> links;
 public:
    void update(BRFoldingTree& tree, const vector<int>& ids) {
        for (int i=0;i+1<tree.seqLen;++i) {
            if (!tree.connectToDownstream[i]) continue;
            BRNode* a=tree.nodes[i]; BRNode* b=tree.nodes[i+1];
            pair<int,int> key(ids[i],ids[i+1]);
            auto found=links.find(key);
            if (found==links.end()) {
                PhophateGroupLocal local=tree.et->pb->getPhoLocal(
                    a->cs2,b->riboAtomCoords,a->rot->improper,b->rot->improper);
                found=links.emplace(key,local).first;
            }
            a->phoLocal=found->second;
            a->pho=PhophateGroup(a->phoLocal,a->cs2);
            a->phoLocalTmp=a->phoLocal;
            a->phoTmp=a->pho;
        }
    }
    size_t size() const { return links.size(); }
};

// Emit molecule-wide sums; the seven identity fields join variant rows to
// contact_variants.csv. The reference row has index -1 and empty identity fields.
static void writeEnergySums(ostream& out, const string& state, long index,
                            const vector<string>& identity, const EnergySums& sums, double total) {
    if (identity.size()!=7 || !sums.agrees(total)) throw runtime_error("Invalid energy component row");
    out<<state<<","<<index;
    for (const string& field : identity) out<<","<<csvQuoted(field);
    out<<setprecision(17)<<","<<sums.bb<<","<<sums.bo<<","<<sums.oo<<","<<sums.clash
       <<","<<sums.pairTotal()<<","<<sums.baseRibose<<","<<sums.basePhosphate
       <<","<<sums.ribosePhosphate<<","<<sums.rotamer<<","<<sums.phosphate
       <<","<<sums.nonpairTotal()<<","<<total<<"\n";
    if (!out) throw runtime_error("Cannot write energy components");
}

static int runMutationTotal(int argc, char** argv) {
    const auto started=chrono::steady_clock::now();
    if (argc<5 || argc>6 || (argc==6 && string(argv[5])!="--verify-full")) {
        cerr<<"Usage: BRiQ_Energy --mutate-total REFERENCE.pdb GEOMETRY.tsv OUTPUT_DIR [--verify-full]"<<endl;
        return 1;
    }
    const bool verify=argc==6;
    const string referencePath=absoluteMutationPath(argv[2]);
    const string geometryPath=absoluteMutationPath(argv[3]);
    if (referencePath.size()<4 || referencePath.substr(referencePath.size()-4)!=".pdb")
        throw runtime_error("Mutation reference must be a .pdb file");
    ifstream raw(referencePath.c_str());
    string rawLine; size_t models=0;
    while (getline(raw,rawLine)) if (rawLine.compare(0,6,"MODEL ")==0) ++models;
    if (!raw.eof() || models>1) throw runtime_error("Select exactly one PDB model");
    MutationGeometry geometry(geometryPath);
    RNAPDB pdb(referencePath,"xxxx");
    RnaAtomLib atomLibrary;
    vector<RNABase*> bases=pdb.getValidBaseList(&atomLibrary);
    if (bases.empty() || bases.size()!=pdb.getBaseList().size())
        throw runtime_error("Reference must contain valid RNA residues");
    validateMutationTopology(bases,geometry);
    AssignRNASS assignment(&pdb,&atomLibrary);
    vector<pair<int,int> > contacts;
    vector<string> contactTypes;
    for (int i=0;i<assignment.len;++i) for (int j=i+1;j<assignment.len;++j) {
        string label=classifyPairContact(&assignment,i,j);
        if (!label.empty() && label!="unannotated") {
            contacts.push_back(make_pair(i,j)); contactTypes.push_back(label);
            for (int index : {i,j}) for (char base : string("ACGU"))
                if (!geometry.bases.count(MutationKey(bases[index]->chainID,bases[index]->baseID,base)))
                    throw runtime_error("Missing one of four variants for a selected contact residue");
        }
    }
    // All lightweight validation precedes output creation and expensive tables.
    if (mkdir(argv[4],0755)!=0 && errno!=EEXIST) throw runtime_error("Cannot create output directory");
    const string output=absoluteMutationPath(argv[4]);
    DIR* directory=opendir(output.c_str());
    if (!directory) throw runtime_error("Cannot inspect output directory");
    bool empty=true; struct dirent* entry;
    while ((entry=readdir(directory))!=NULL)
        if (string(entry->d_name)!="." && string(entry->d_name)!="..") empty=false;
    closedir(directory);
    if (!empty) throw runtime_error("Mutation output directory must be empty");
    const string fingerprint=batchFingerprint({referencePath,geometryPath});
    const auto loadStart=chrono::steady_clock::now();
    cerr<<"[mutate-total] Loading shared scoring libraries once"<<endl;
    EnergyTable tables;
    BaseRotamerLib rotamers;
    const double librarySeconds=chrono::duration<double>(chrono::steady_clock::now()-loadStart).count();
    unique_ptr<BRFoldingTree> tree;
    double nativeEnergy;
    {
        MutationWorkingDirectory cwd;
        vector<string> args={"BRiQ_Energy",referencePath,"--output-dir",output,
            "--seed","0","--weights","1","1","1","--pair-ownership","pdb"};
        vector<char*> pointers;
        for (string& arg : args) pointers.push_back(&arg[0]);
        if (runEnergy(args.size(),pointers.data(),&tables,&nativeEnergy,&rotamers,
                NULL,NULL,true,NULL,&tree)!=0)
            throw runtime_error("Reference scoring failed");
    }
    if (tree->seqLen!=static_cast<int>(bases.size())) throw runtime_error("Reference length mismatch");
    IncrementalEnergy incremental;
    EnergySums nativeSums;
    const double checkedNative=incremental.score(*tree,&nativeSums);
    if (fabs(checkedNative-nativeEnergy)>1e-9) throw runtime_error("Native component scorer mismatch");
    const auto cacheStart=chrono::steady_clock::now();
    RNABaseLib standardizer;
    map<MutationKey,CachedMutationNode> nodes;
    map<MutationKey,int> nodeIds;
    for (auto& entry : geometry.bases) {
        size_t index=0;
        while (bases[index]->chainID!=get<0>(entry.first) || bases[index]->baseID!=get<1>(entry.first)) ++index;
        unique_ptr<RNABase> standard(standardizer.toStandardBase(entry.second.get()));
        standard->baseSeqID=index;
        CachedMutationNode cached;
        cached.rotamer.reset(standard->backboneComplete() ? new BaseRotamer(standard.get()) :
            new BaseRotamer(*rotamers.getLowestEnergyRotamer(standard->baseTypeInt)));
        cached.node.reset(new BRNode(standard.get(),cached.rotamer.get()));
        cached.node->connectToNeighbor=tree->connectToDownstream[index];
        nodeIds[entry.first]=tree->seqLen+nodes.size();
        nodes.emplace(entry.first,move(cached));
    }
    const double cacheSeconds=chrono::duration<double>(chrono::steady_clock::now()-cacheStart).count();
    MutationNodeGuard guard(*tree);
    MutationPhosphates phosphates;
    vector<int> ids(tree->seqLen);
    const string temporary=output+"/contact_variants.csv.tmp";
    ofstream totals(temporary.c_str());
    if (!totals) throw runtime_error("Cannot create variant totals");
    totals<<setprecision(17)<<"chain_i,residue_i,chain_j,residue_j,contact_type,base_i,base_j,native_E_total,E_total,delta_E_total\n";
    const string componentsTemporary=output+"/energy_components.csv.tmp";
    ofstream components(componentsTemporary.c_str());
    if (!components) throw runtime_error("Cannot create energy component summary");
    components<<"state,variant_index,chain_i,residue_i,chain_j,residue_j,contact_type,base_i,base_j,"
              "E_bb,E_bo,E_oo,E_clash,E_pair,E_base_ribose,E_base_phosphate,E_ribose_phosphate,"
              "E_rotamer,E_phosphate_local,E_nonpair,E_total\n";
    writeEnergySums(components,"reference",-1,vector<string>(7,""),nativeSums,nativeEnergy);
    const auto variantsStart=chrono::steady_clock::now();
    size_t count=0;
    for (size_t contact=0;contact<contacts.size();++contact) {
        const int i=contacts[contact].first, j=contacts[contact].second;
        cerr<<"[mutate-total "<<contact+1<<"/"<<contacts.size()<<"] "
            <<bases[i]->chainID<<":"<<bases[i]->baseID<<" "<<bases[j]->chainID<<":"<<bases[j]->baseID<<endl;
        for (char a : string("ACGU")) for (char b : string("ACGU")) {
            guard.restore();
            for (int k=0;k<tree->seqLen;++k) ids[k]=k;
            MutationKey first(bases[i]->chainID,bases[i]->baseID,a);
            MutationKey second(bases[j]->chainID,bases[j]->baseID,b);
            tree->nodes[i]=nodes.at(first).node.get(); ids[i]=nodeIds.at(first);
            tree->nodes[j]=nodes.at(second).node.get(); ids[j]=nodeIds.at(second);
            tree->seq[i]=tree->nodes[i]->baseType; tree->seq[j]=tree->nodes[j]->baseType;
            phosphates.update(*tree,ids);
            EnergySums sums;
            const double energy=incremental.score(*tree,&sums);
            if (!isfinite(energy)) throw runtime_error("Non-finite mutant energy");
            if (verify) {
                // Independently rebuild every phosphate and evaluate all terms.
                tree->updatePhoGroups();
                std::array<double,9> fullSums;
                double full=tree->totalEnergy(1,1,1,0,false,NULL,NULL,
                    vector<PairEnergyMetadata>(),vector<string>(),true,&fullSums);
                const auto cachedSums=sums.values();
                for (size_t component=0;component<fullSums.size();++component)
                    if (!isfinite(fullSums[component]) || fabs(fullSums[component]-cachedSums[component])>1e-9)
                        throw runtime_error("Full/memory component mismatch at variant " + to_string(count));
                if (!isfinite(full) || fabs(full-energy)>1e-9)
                    throw runtime_error("Full/memory score mismatch at variant " + to_string(count));
            }
            totals<<csvQuoted(bases[i]->chainID)<<","<<csvQuoted(bases[i]->baseID)<<","
                <<csvQuoted(bases[j]->chainID)<<","<<csvQuoted(bases[j]->baseID)<<","
                <<csvQuoted(contactTypes[contact])<<","<<a<<","<<b<<","<<nativeEnergy<<","
                <<energy<<","<<energy-nativeEnergy<<"\n";
            if (!totals) throw runtime_error("Cannot write mutant energy");
            writeEnergySums(components,"variant",count,
                {bases[i]->chainID,bases[i]->baseID,bases[j]->chainID,bases[j]->baseID,
                 contactTypes[contact],string(1,a),string(1,b)},sums,energy);
            ++count;
        }
    }
    guard.restore();
    totals.close();
    if (!totals.good() || rename(temporary.c_str(),(output+"/contact_variants.csv").c_str())!=0)
        throw runtime_error("Cannot finalize mutant totals");
    components.close();
    if (!components.good() || rename(componentsTemporary.c_str(),(output+"/energy_components.csv").c_str())!=0)
        throw runtime_error("Cannot finalize energy components");
    const double variantsSeconds=chrono::duration<double>(chrono::steady_clock::now()-variantsStart).count();
    ofstream report((output+"/mutation_run.json").c_str());
    report<<setprecision(17)<<"{\n  \"source_version\": \""<<BRIQ_SOURCE_VERSION<<"\",\n"
        <<"  \"engine\": \"memory\",\n  \"geometry_protocol\": \"supplied_coordinates\",\n"
        <<"  \"seed\": 0,\n  \"weights\": [1, 1, 1],\n"
        <<"  \"input_fingerprint\": \""<<fingerprint<<"\",\n"
        <<"  \"energy_components_file\": \"energy_components.csv\",\n"
        <<"  \"contacts\": "<<contacts.size()<<",\n  \"variants\": "<<count<<",\n"
        <<"  \"native_energy\": "<<nativeEnergy<<",\n"
        <<"  \"verify_full\": "<<(verify ? "true" : "false")<<",\n"
        <<"  \"phosphate_cache_entries\": "<<phosphates.size()<<",\n"
        <<"  \"library_seconds\": "<<librarySeconds<<",\n"
        <<"  \"node_cache_seconds\": "<<cacheSeconds<<",\n"
        <<"  \"variant_seconds\": "<<variantsSeconds<<",\n"
        <<"  \"total_seconds\": "<<chrono::duration<double>(chrono::steady_clock::now()-started).count()<<"\n}\n";
    report.close();
    if (!report.good()) throw runtime_error("Cannot write mutation run metadata");
    cerr<<"[mutate-total] Complete: "<<count<<" variants"<<endl;
    return 0;
}
} // namespace
#endif
