/*
 * RunEnergy.cpp
 *
 *  Created on: 2025-04-28
 *      Author: tzok
 */

#include "pred/BRFoldingTree.h"
#include "pred/IncrementalEnergy.h"
#include "model/AssignRNASS.h"
#include "model/ProteinRep.h"
#include <fstream>
#include <cstdio>
#include <fcntl.h>
#include <stdexcept>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <memory>
#include <utility>
#include <chrono>
#include <sstream>
#include <set>
#include <vector>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>

#ifndef BRIQ_SOURCE_VERSION
#define BRIQ_SOURCE_VERSION "unknown"
#endif

using namespace NSPpred;
using namespace std;

static string trimAtomName(const string &name) {
  size_t first = name.find_first_not_of(" \t");
  if (first == string::npos)
    return "";
  size_t last = name.find_last_not_of(" \t");
  return name.substr(first, last - first + 1);
}

static bool isHeavyAtom(const Atom* atom) {
  string element = trimAtomName(atom->type);
  string atomName = trimAtomName(atom->name);
  if (!element.empty() &&
      (toupper(element[0]) == 'H' || toupper(element[0]) == 'D'))
    return false;
  if (element.empty() && !atomName.empty() &&
      (atomName[0] == 'H' || atomName[0] == 'D'))
    return false;
  return true;
}

static double minimumAtomDistance(const vector<Atom *> &first,
                                  const vector<Atom *> &second,
                                  const vector<string> *allowedNames = NULL) {
  double minimum = numeric_limits<double>::max();
  for (size_t i = 0; i < first.size(); i++) {
    Atom *a = first[i];
    if (!isHeavyAtom(a))
      continue;
    for (size_t j = 0; j < second.size(); j++) {
      Atom *b = second[j];
      if (!isHeavyAtom(b))
        continue;
      if (allowedNames != NULL) {
        string name = trimAtomName(b->name);
        if (find(allowedNames->begin(), allowedNames->end(), name) ==
            allowedNames->end())
          continue;
      }
      double distance = a->coord.distance(b->coord);
      if (distance < minimum)
        minimum = distance;
    }
  }
  return minimum;
}

static bool isStackingContact(RNABase *first, RNABase *second) {
  XYZ firstNormal = first->getBaseNormVector();
  XYZ secondNormal = second->getBaseNormVector();
  if (firstNormal.length() == 0.0 || secondNormal.length() == 0.0)
    return false;

  double normalDot = firstNormal.x_ * secondNormal.x_ +
                     firstNormal.y_ * secondNormal.y_ +
                     firstNormal.z_ * secondNormal.z_;
  double centerDistance = first->distanceTo(*second);
  return fabs(normalDot) >= cos(30.0 * 3.14159265358979323846 / 180.0) &&
         first->planeDistance(*second) <= 4.0 && centerDistance <= 8.0;
}

static string classifyPairContact(AssignRNASS *assignment, int i, int j) {
  RNABase *first = assignment->baseList[i];
  RNABase *second = assignment->baseList[j];
  vector<string> labels;
  if (assignment->pairIndex[i] == j || assignment->pairIndex[j] == i)
    labels.push_back("base_pair");
  if (assignment->nwcPairIndex[i] == j ||
      assignment->nwcPairIndex[j] == i)
    labels.push_back("non_watson_crick_pair");

  bool neighbors = j == i + 1 && first->chainID == second->chainID &&
                   first->connectToNeighbor(*second);
  double baseDistance = minimumAtomDistance(*first->getSidechainAtoms(),
                                             *second->getSidechainAtoms());
  if (neighbors && isStackingContact(first, second))
    labels.push_back("stacking");
  else if (baseDistance <= 4.5)
    labels.push_back("base_contact");

  const vector<string> sugarNames = {"C1'", "C2'", "C3'", "C4'", "C5'",
                                     "O2'", "O3'", "O4'", "O5'"};
  const vector<string> phosphateNames = {"P", "OP1", "OP2", "OP3",
                                          "O1P", "O2P"};
  double baseSugarDistance = minimumAtomDistance(
      *first->getSidechainAtoms(), *second->getBackboneAtoms(), &sugarNames);
  baseSugarDistance = min(
      baseSugarDistance,
      minimumAtomDistance(*second->getSidechainAtoms(),
                          *first->getBackboneAtoms(), &sugarNames));
  if (baseSugarDistance <= 4.0)
    labels.push_back("base_ribose");

  double basePhosphateDistance = minimumAtomDistance(
      *first->getSidechainAtoms(), *second->getBackboneAtoms(),
      &phosphateNames);
  basePhosphateDistance = min(
      basePhosphateDistance,
      minimumAtomDistance(*second->getSidechainAtoms(),
                          *first->getBackboneAtoms(), &phosphateNames));
  if (basePhosphateDistance <= 4.0)
    labels.push_back("base_phosphate");

  if (labels.empty())
    return "unannotated";
  ostringstream result;
  for (size_t k = 0; k < labels.size(); k++) {
    if (k > 0)
      result << ";";
    result << labels[k];
  }
  return result.str();
}

static bool parseDouble(const char *text, double &value) {
  char *end = NULL;
  value = strtod(text, &end);
  return end != text && *end == '\0' && isfinite(value);
}

#include "tools/ScoringCheckpoint.h"

bool checkInputFile(const string &inputFile) {
  ifstream file;
  cout << "check input file: " << endl;
  file.open(inputFile.c_str(), ios::in);
  if (!file.is_open()) {
    cout << "fail to open inputFile: " << inputFile << endl;
    return false;
  }
  file.close();

  NSPtools::InputParser input(inputFile);
  string pdbFile = input.getValue("pdb");
  string baseSeq = input.getValue("seq");
  string baseSec = input.getValue("sec");
  string nwcSec = input.getValue("nwc");

  if (pdbFile == "") {
    cout << "invalid inputFile: initial pdb required" << endl;
    return false;
  }

  if (baseSeq == "") {
    cout << "invalid inputFile: base sequence required" << endl;
    return false;
  }

  if (baseSec == "") {
    cout << "invalid inputFile: Watson-Crick pairing sequence required" << endl;
    return false;
  }

  if (nwcSec == "") {
    cout << "invalid inputFile: non-Watson-Crick pairing sequence required"
         << endl;
    return false;
  }

  if (baseSec.length() != baseSeq.length()) {
    cout << "the length of Watson-Crick pairing sequence must equal to the "
            "length of base sequence"
         << endl;
    return false;
  }

  if (nwcSec.length() != baseSeq.length()) {
    cout << "the length of non-Watson-Crick pairing sequence must equal to the "
            "length of base sequence"
         << endl;
    return false;
  }

  file.open(pdbFile.c_str(), ios::in);
  if (!file.is_open()) {
    cout << "fail to open pdb file: " << pdbFile << endl;
    return false;
  }
  return true;
}

void printHelp() {
  cout << "       BRiQ_Energy --mutate-total REFERENCE.pdb GEOMETRY.tsv OUTPUT_DIR [--verify-full]" << endl;
  cout << "       BRiQ_Energy --batch-total JOBS.tsv TOTALS.csv (one absolute PDB path per row)" << endl;
  cout << "         [--engine legacy|shared|score|incremental] [--profile]"
          " [--resume|--resume-verify] [--flush-every N]" << endl;
  cout << "       BRiQ_Energy --batch JOBS.tsv (seed 0, weights 1/1/1, PDB ownership)" << endl;
  cout << "Usage: BRiQ_Energy INPUT_FILE [--pair-energies PAIRS.csv]"
          " [--nonpair NONPAIR.csv] [--components COMPONENTS.tsv]"
          " [--weights BREAK CONNECT CLASH] [--pair-ownership pdb|node]"
          " [--seed SEED]" << endl;
  cout << "       BRiQ_Energy --prepare INPUT.pdb OUTPUT_FILE" << endl;
  cout << "       BRiQ_Energy INPUT.pdb --output-dir DIRECTORY"
          " [--weights BREAK CONNECT CLASH] [--pair-ownership pdb|node]"
          " [--seed SEED]" << endl;
}

bool prepareEnergyInput(const string &pdbFile, const string &outputFile) {
  if (pdbFile.empty() || outputFile.empty()) {
    cerr << "PDB input and energy input output paths must not be empty" << endl;
    return false;
  }

  ifstream pdbCheck(pdbFile.c_str(), ios::in);
  if (!pdbCheck.is_open()) {
    cerr << "fail to open pdb file: " << pdbFile << endl;
    return false;
  }
  pdbCheck.close();

  RNAPDB *pdb = new RNAPDB(pdbFile, "xxxx");
  RnaAtomLib *atomLibrary = new RnaAtomLib();
  AssignRNASS *assignment = new AssignRNASS(pdb, atomLibrary);

  ofstream output(outputFile.c_str(), ios::out);
  if (!output.is_open()) {
    cerr << "fail to open output file: " << outputFile << endl;
    delete assignment;
    delete atomLibrary;
    delete pdb;
    return false;
  }

  output << "pdb " << pdbFile << endl;
  output << "seq " << assignment->seq << endl;
  output << "sec " << assignment->ssSeq << endl;
  output << "nwc " << assignment->nwcSeq << endl;
  set<int> breaks(assignment->breakList.begin(), assignment->breakList.end());
  for (int i = 0; i + 1 < assignment->len; i++) {
    if (assignment->baseList[i]->chainID != assignment->baseList[i + 1]->chainID)
      breaks.insert(i);
  }
  if (!breaks.empty()) {
    output << "break";
    for (set<int>::const_iterator it = breaks.begin(); it != breaks.end(); ++it)
      output << " " << *it;
    output << endl;
  }
  output.close();

  delete assignment;
  delete atomLibrary;
  delete pdb;
  return output.good();
}

// For fixed equal-weight scoring, secondary-structure/sampling annotations do
// not affect the score. Preserve exactly the same valid residues and breaks.
static bool prepareScoringInput(const string& pdbPath, const string& outputPath) {
    RNAPDB pdb(pdbPath, "xxxx");
    RnaAtomLib atomLibrary;
    vector<RNABase*> bases=pdb.getValidBaseList(&atomLibrary);
    if (bases.empty() || bases.size()!=pdb.getBaseList().size()) {
        cerr << "Fast scoring requires a nonempty, valid RNA input" << endl;
        return false;
    }
    string sequence;
    for (size_t i=0;i<bases.size();++i) sequence+=bases[i]->baseType;
    ofstream output(outputPath.c_str());
    output << "pdb " << pdbPath << "\nseq " << sequence
           << "\nsec " << string(bases.size(),'.')
           << "\nnwc " << string(bases.size(),'.') << "\n";
    vector<size_t> breaks;
    for (size_t i=0;i+1<bases.size();++i)
        if (bases[i]->chainID!=bases[i+1]->chainID ||
            !bases[i]->connectToNeighbor(*bases[i+1])) breaks.push_back(i);
    if (!breaks.empty()) {
        output << "break";
        for (size_t i=0;i<breaks.size();++i) output << " " << breaks[i];
        output << "\n";
    }
    output.close();
    return output.good();
}

static bool startOutputDirectory(const string &source, const string &directory,
                                 string &preparedFile) {
  size_t separator = source.find_last_of("/\\");
  string pdbName = source.substr(separator == string::npos ? 0 : separator + 1);
  size_t extension = pdbName.find_last_of('.');
  string stem = pdbName.substr(0, extension);
  preparedFile = stem + ".briq";
  if (stem.empty() || pdbName == "pairs.csv" || pdbName == "nonpair.csv" ||
      pdbName == "run.json" || pdbName == preparedFile) {
    cerr << "input filename conflicts with output filenames" << endl;
    return false;
  }
  ifstream input(source.c_str(), ios::binary);
  if (!input) {
    cerr << "cannot open input PDB: " << source << endl;
    return false;
  }
  if (mkdir(directory.c_str(), 0755) != 0 && errno != EEXIST) {
    cerr << "cannot create output directory (parent must exist): " << directory << endl;
    return false;
  }
  DIR *dir = opendir(directory.c_str());
  if (!dir) {
    cerr << "cannot open output directory: " << directory << endl;
    return false;
  }
  bool empty = true;
  struct dirent *entry;
  while ((entry = readdir(dir)) != NULL) {
    string name(entry->d_name);
    if (name != "." && name != "..") empty = false;
  }
  closedir(dir);
  if (!empty) {
    cerr << "output directory must be empty: " << directory << endl;
    return false;
  }
  if (chdir(directory.c_str()) != 0) {
    cerr << "cannot enter output directory" << endl;
    return false;
  }
  ofstream copy(pdbName.c_str(), ios::binary);
  copy << input.rdbuf();
  copy.close();
  if (!copy.good()) {
    cerr << "cannot copy input PDB" << endl;
    return false;
  }
  return prepareEnergyInput(pdbName, preparedFile);
}

// Fast scoring is serial. Buffer verbose upstream diagnostics and display them
// only on failure; keep progress/profile lines on stderr. This avoids thousands
// of flushed node messages for every variant without hiding error diagnostics.
class ScoringDiagnostics {
    ostringstream captured;
    streambuf* previous;
    bool success;
 public:
    explicit ScoringDiagnostics(bool quiet): previous(NULL), success(false) {
        if (quiet) previous=cout.rdbuf(captured.rdbuf());
    }
    void complete() {success=true;}
    ~ScoringDiagnostics() {
        if (previous != NULL) {
            cout.rdbuf(previous);
            if (!success) cerr<<captured.str();
        }
    }
};

static int runEnergy(int argc, char **argv, EnergyTable *sharedTables = NULL,
                     double *totalResult = NULL, BaseRotamerLib *rotamers = NULL,
                     FragmentLibrary *fragments = NULL, ScoringProfile *profile = NULL,
                     bool scoringOnly = false, IncrementalEnergy *incremental = NULL,
                     unique_ptr<BRFoldingTree> *retainedTree = NULL) {
  if (argc > 1 && string(argv[1]) == "--prepare") {
    if (argc != 4) {
      printHelp();
      return 1;
    }
    return prepareEnergyInput(argv[2], argv[3]) ? 0 : 1;
  }

  if (argc == 2 && (string(argv[1]) == "-h" ||
                    string(argv[1]) == "--help")) {
    printHelp();
    return 0;
  }
  if (argc < 2) {
    printHelp();
    return 1;
  }

  string inputFile = string(argv[1]);
  string outputDirectory;
  string pairFile;
  string nonpairFile;
  string componentFile;
  double breakWeight = 1.0;
  double connectWeight = 1.0;
  double clashWeight = 1.0;
  bool pdbOwnership = true;
  unsigned int randomSeed = static_cast<unsigned int>(time(NULL));
  for (int arg = 2; arg < argc; arg++) {
    string option(argv[arg]);
    if (option == "--output-dir") {
      if (arg + 1 >= argc || !outputDirectory.empty() || !argv[arg + 1][0]) {
        cerr << "--output-dir requires one nonempty directory" << endl;
        return 1;
      }
      outputDirectory = argv[++arg];
    } else if (option == "--pair-energies") {
      if (arg + 1 >= argc || !pairFile.empty()) {
        cerr << "missing or duplicate pair-energy output path" << endl;
        return 1;
      }
      pairFile = argv[++arg];
      if (pairFile.empty()) {
        cerr << "pair-energy output path must not be empty" << endl;
        return 1;
      }
    } else if (option == "--components") {
      if (arg + 1 >= argc || !componentFile.empty()) {
        cerr << "missing or duplicate component output path" << endl;
        return 1;
      }
      componentFile = argv[++arg];
      if (componentFile.empty()) {
        cerr << "component output path must not be empty" << endl;
        return 1;
      }
    } else if (option == "--nonpair") {
      if (arg + 1 >= argc || !nonpairFile.empty()) {
        cerr << "missing or duplicate nonpair output path" << endl;
        return 1;
      }
      nonpairFile = argv[++arg];
      if (nonpairFile.empty()) {
        cerr << "nonpair output path must not be empty" << endl;
        return 1;
      }
    } else if (option == "--pair-ownership") {
      if (arg + 1 >= argc) {
        cerr << "--pair-ownership requires pdb or node" << endl;
        return 1;
      }
      string ownership(argv[++arg]);
      if (ownership != "pdb" && ownership != "node") {
        cerr << "--pair-ownership requires pdb or node" << endl;
        return 1;
      }
      pdbOwnership = ownership == "pdb";
    } else if (option == "--seed") {
      if (arg + 1 >= argc) {
        cerr << "--seed requires an unsigned integer" << endl;
        return 1;
      }
      string seedText(argv[++arg]);
      if (seedText.empty() || seedText.find_first_not_of("0123456789") != string::npos) {
        cerr << "--seed requires an unsigned integer" << endl;
        return 1;
      }
      char *end = NULL;
      errno = 0;
      unsigned long parsedSeed = strtoul(seedText.c_str(), &end, 10);
      if (errno == ERANGE || *end != '\0' || parsedSeed > numeric_limits<unsigned int>::max()) {
        cerr << "--seed is outside the unsigned integer range" << endl;
        return 1;
      }
      randomSeed = static_cast<unsigned int>(parsedSeed);
    } else if (option == "--weights") {
      if (arg + 3 >= argc ||
          !parseDouble(argv[arg + 1], breakWeight) ||
          !parseDouble(argv[arg + 2], connectWeight) ||
          !parseDouble(argv[arg + 3], clashWeight)) {
        cerr << "--weights requires three numeric values: BREAK CONNECT CLASH"
             << endl;
        return 1;
      }
      arg += 3;
    } else {
      cerr << "unknown option: " << option << endl;
      printHelp();
      return 1;
    }
  }
  if (!outputDirectory.empty()) {
    if (!pairFile.empty() || !nonpairFile.empty() || !componentFile.empty()) {
      cerr << "--output-dir cannot be combined with individual output paths" << endl;
      return 1;
    }
    string preparedFile;
    if (!startOutputDirectory(inputFile, outputDirectory, preparedFile)) return 1;
    inputFile = preparedFile;
    pairFile = "pairs.csv";
    nonpairFile = "nonpair.csv";
  }
  bool exportPairs = !pairFile.empty();
  if (!exportPairs && !nonpairFile.empty()) {
    cerr << "--nonpair requires --pair-energies" << endl;
    return 1;
  }
  if (exportPairs && nonpairFile.empty())
    nonpairFile = pairFile + ".nonpair.csv";
  if (exportPairs && pairFile == nonpairFile) {
    cerr << "pair and nonpair output paths must be different" << endl;
    return 1;
  }
  if ((!pairFile.empty() && pairFile == componentFile) ||
      (!nonpairFile.empty() && nonpairFile == componentFile)) {
    cerr << "pair, nonpair, and component output paths must be different"
         << endl;
    return 1;
  }

  srand(randomSeed);

  ScoringDiagnostics diagnostics(scoringOnly);
  if (!checkInputFile(inputFile)) return 1;

  unique_ptr<EnergyTable> ownedTables;
  if (sharedTables == NULL) ownedTables.reset(new EnergyTable());
  EnergyTable *et = sharedTables != NULL ? sharedTables : ownedTables.get();
  unique_ptr<BRFoldingTree> ownedTree(new BRFoldingTree(inputFile, et, rotamers, fragments, profile, scoringOnly));
  BRFoldingTree *ft = ownedTree.get();
  const auto scoreStart = chrono::steady_clock::now();
  double energy;
  ofstream components;
  if (!componentFile.empty()) {
    components.open(componentFile.c_str(), ios::out);
    if (!components.is_open()) {
      cerr << "fail to open component output file: " << componentFile << endl;
      return 1;
    }
  }
  if (exportPairs) {
    ofstream pairOutput(pairFile.c_str(), ios::out);
    ofstream nonpairOutput(nonpairFile.c_str(), ios::out);
    if (!pairOutput.is_open() || !nonpairOutput.is_open()) {
      cerr << "fail to open pair or nonpair output file" << endl;
      if (pairOutput.is_open())
        pairOutput.close();
      if (nonpairOutput.is_open())
        nonpairOutput.close();
      return 1;
    }
    pairOutput << setprecision(17);
    nonpairOutput << setprecision(17);

    NSPtools::InputParser input(inputFile);
    string pdbFile = input.getValue("pdb");
    RNAPDB *pdb = new RNAPDB(pdbFile, "xxxx");
    RnaAtomLib *atomLibrary = new RnaAtomLib();
    AssignRNASS *assignment = new AssignRNASS(pdb, atomLibrary);
    if (assignment->len != ft->seqLen) {
      cerr << "PDB residue count does not match the BRiQ input sequence length"
           << endl;
      delete assignment;
      delete atomLibrary;
      delete pdb;
      pairOutput.close();
      nonpairOutput.close();
      return 1;
    }
    vector<PairEnergyMetadata> residueMetadata;
    vector<string> contactTypes((size_t)assignment->len * assignment->len);
    for (int i = 0; i < assignment->len; i++) {
      PairEnergyMetadata metadata;
      metadata.chain = assignment->baseList[i]->chainID;
      metadata.residue = assignment->baseList[i]->baseID;
      residueMetadata.push_back(metadata);
      for (int j = i + 1; j < assignment->len; j++) {
        string contact = classifyPairContact(assignment, i, j);
        contactTypes[(size_t)i * assignment->len + j] = contact;
        contactTypes[(size_t)j * assignment->len + i] = contact;
      }
    }

    double exportedEnergy = ft->totalEnergy(
        breakWeight, connectWeight, clashWeight, 0.0, false, &pairOutput,
        &nonpairOutput, residueMetadata, contactTypes, pdbOwnership);
    if (components.is_open()) {
      energy = ft->totalEnergy(breakWeight, connectWeight, clashWeight, 0.0,
                               false, &components);
      if (fabs(energy - exportedEnergy) > 1e-9) {
        cerr << "pair-energy export does not reproduce original totalEnergy: "
             << setprecision(17) << exportedEnergy << " vs " << energy << endl;
        return 1;
      }
    } else {
      energy = exportedEnergy;
    }
    pairOutput.close();
    nonpairOutput.close();
    delete assignment;
    delete atomLibrary;
    delete pdb;
    if (!pairOutput.good() || !nonpairOutput.good()) {
      cerr << "failed while writing pair-energy output files" << endl;
      return 1;
    }
  } else if (components.is_open()) {
    energy = ft->totalEnergy(breakWeight, connectWeight, clashWeight, 0.0,
                             false, &components);
  } else if (incremental != NULL) {
    energy = incremental->score(*ft);
  } else {
    energy = ft->totalEnergy(breakWeight, connectWeight, clashWeight, 0.0,
                             false);
  }
  if (components.is_open()) {
    components.close();
    if (!components.good()) {
      cerr << "failed while writing component output file: " << componentFile
           << endl;
      return 1;
    }
  }
  if (!isfinite(energy)) {
    cerr << "BRiQ produced a non-finite energy" << endl;
    return 1;
  }
  if (profile != NULL) profile->scoreSeconds = chrono::duration<double>(chrono::steady_clock::now()-scoreStart).count();
  if (totalResult != NULL) *totalResult = energy;
  if (!outputDirectory.empty()) {
    ofstream run("run.json");
    run << setprecision(17)
        << "{\n  \"program\": \"BRiQ_Energy\",\n"
        << "  \"version\": \"1.0\",\n"
        << "  \"exporter_version\": \"0.3.0\",\n"
        << "  \"source_version\": \"" << BRIQ_SOURCE_VERSION << "\",\n"
        << "  \"parameters\": {\n"
        << "    \"break_weight\": " << breakWeight << ",\n"
        << "    \"connect_weight\": " << connectWeight << ",\n"
        << "    \"clash_weight\": " << clashWeight << ",\n"
        << "    \"seed\": " << randomSeed << ",\n"
        << "    \"pair_ownership\": \"" << (pdbOwnership ? "pdb" : "node") << "\"\n"
        << "  },\n  \"energy_total\": " << energy << "\n}\n";
    run.close();
    if (!run.good()) {
      cerr << "cannot write run.json" << endl;
      return 1;
    }
  }
  if (exportPairs)
    cerr << "pair_ownership " << (pdbOwnership ? "pdb" : "node") << endl;
  if (!scoringOnly) cerr << inputFile << " " << setprecision(17) << energy << endl;
  if (retainedTree != NULL) *retainedTree = move(ownedTree);
  diagnostics.complete();
  return 0;
}

struct BatchJob {
  string input;
  string output;
};

// Validate every job before loading tables or creating any job output.
static bool readBatchJobs(const string &manifestPath, vector<BatchJob> &jobs) {
  ifstream manifest(manifestPath.c_str());
  if (!manifest) {
    cerr << "Cannot open batch manifest: " << manifestPath << endl;
    return false;
  }
  set<string> destinations;
  string line;
  size_t lineNumber = 0;
  while (getline(manifest, line)) {
    ++lineNumber;
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty()) continue;
    size_t separator = line.find('\t');
    if (separator == string::npos ||
        line.find('\t', separator + 1) != string::npos) {
      cerr << "Manifest line " << lineNumber
           << ": expected PDB<TAB>OUTPUT_DIRECTORY (no header)" << endl;
      return false;
    }
    BatchJob job = {line.substr(0, separator), line.substr(separator + 1)};
    if (job.input.empty() || job.output.empty() ||
        job.input[0] != '/' || job.output[0] != '/') {
      cerr << "Manifest line " << lineNumber
           << ": paths must be absolute container paths" << endl;
      return false;
    }
    struct stat info;
    if (stat(job.input.c_str(), &info) != 0 || !S_ISREG(info.st_mode)) {
      cerr << "Cannot read batch PDB: " << job.input << endl;
      return false;
    }
    ifstream input(job.input.c_str());
    if (!input) {
      cerr << "Cannot open batch PDB: " << job.input << endl;
      return false;
    }
    while (job.output.size() > 1 && job.output.back() == '/') job.output.pop_back();
    size_t slash = job.output.find_last_of('/');
    string parent = slash == 0 ? "/" : job.output.substr(0, slash);
    string name = job.output.substr(slash + 1);
    if (name.empty() || name == "." || name == ".." ||
        stat(parent.c_str(), &info) != 0 || !S_ISDIR(info.st_mode)) {
      cerr << "Batch output parent must exist: " << job.output << endl;
      return false;
    }
    char *resolvedParent = realpath(parent.c_str(), NULL);
    if (resolvedParent == NULL) {
      cerr << "Cannot resolve batch output parent: " << parent << endl;
      return false;
    }
    job.output = string(resolvedParent) + "/" + name;
    free(resolvedParent);
    if (stat(job.output.c_str(), &info) == 0) {
      if (!S_ISDIR(info.st_mode)) {
        cerr << "Batch output is not a directory: " << job.output << endl;
        return false;
      }
      char *resolved = realpath(job.output.c_str(), NULL);
      if (resolved == NULL) return false;
      job.output = resolved;
      free(resolved);
      DIR *directory = opendir(job.output.c_str());
      if (directory == NULL) return false;
      bool empty = true;
      struct dirent *entry;
      while ((entry = readdir(directory)) != NULL) {
        string entryName(entry->d_name);
        if (entryName != "." && entryName != "..") empty = false;
      }
      closedir(directory);
      if (!empty) {
        cerr << "Batch output directory must be empty: " << job.output << endl;
        return false;
      }
    } else if (errno != ENOENT) {
      cerr << "Cannot access batch output directory: " << job.output << endl;
      return false;
    }
    if (!destinations.insert(job.output).second) {
      cerr << "Duplicate batch output directory: " << job.output << endl;
      return false;
    }
    jobs.push_back(job);
  }
  if (!manifest.eof()) {
    cerr << "Failed while reading batch manifest" << endl;
    return false;
  }
  if (jobs.empty()) {
    cerr << "Batch manifest is empty" << endl;
    return false;
  }
  return true;
}

static bool completeBatchOutput(const BatchJob &job) {
  size_t separator = job.input.find_last_of('/');
  string name = job.input.substr(separator + 1);
  string stem = name.substr(0, name.find_last_of('.'));
  const vector<string> expected = {
    name, stem + ".briq", "pairs.csv", "nonpair.csv", "run.json"
  };
  for (size_t i = 0; i < expected.size(); ++i) {
    ifstream file((job.output + "/" + expected[i]).c_str());
    if (!file) {
      cerr << "Missing batch output: " << job.output << "/" << expected[i] << endl;
      return false;
    }
  }
  return true;
}

static int runBatch(int argc, char **argv) {
  if (argc != 3) {
    cerr << "Usage: BRiQ_Energy --batch JOBS.tsv" << endl;
    return 1;
  }
  vector<BatchJob> jobs;
  if (!readBatchJobs(argv[2], jobs)) return 1;

  char *directory = realpath(".", NULL);
  if (directory == NULL) {
    cerr << "Cannot resolve working directory" << endl;
    return 1;
  }
  string initialDirectory(directory);
  free(directory);

  srand(0);
  const chrono::steady_clock::time_point start = chrono::steady_clock::now();
  cerr << "[batch] Loading shared energy tables" << endl;
  EnergyTable tables;
  cerr << "[batch] Tables ready in "
       << chrono::duration<double>(chrono::steady_clock::now() - start).count()
       << " s" << endl;

  for (size_t index = 0; index < jobs.size(); ++index) {
    cerr << "[batch " << index + 1 << "/" << jobs.size() << "] "
         << jobs[index].input << endl;
    vector<string> arguments = {
      "BRiQ_Energy", jobs[index].input,
      "--output-dir", jobs[index].output,
      "--seed", "0", "--weights", "1", "1", "1",
      "--pair-ownership", "pdb"
    };
    vector<char *> pointers;
    for (size_t i = 0; i < arguments.size(); ++i) pointers.push_back(&arguments[i][0]);
    pointers.push_back(NULL);
    int status = runEnergy(static_cast<int>(arguments.size()), pointers.data(), &tables);
    if (chdir(initialDirectory.c_str()) != 0) {
      cerr << "Cannot restore batch working directory" << endl;
      return 1;
    }
    if (status != 0 || !completeBatchOutput(jobs[index])) {
      cerr << "Batch stopped at job " << index + 1 << endl;
      return status != 0 ? status : 1;
    }
  }
  return 0;
}

// Internal preparation files live in the container, never in the mounted workspace.
// A private symlink avoids copying the PDB or changing the working directory.
class TotalScratch {
 public:
  string directory;
  string pdb;
  string prepared;
  TotalScratch() {
    char pattern[] = "/tmp/briq-total-XXXXXX";
    char *created = mkdtemp(pattern);
    if (created == NULL) throw runtime_error("Cannot create total-scoring scratch directory");
    directory = created;
    pdb = directory + "/input.pdb";
    prepared = directory + "/input.briq";
  }
  ~TotalScratch() {
    unlink(prepared.c_str());
    unlink(pdb.c_str());
    rmdir(directory.c_str());
  }
 private:
  TotalScratch(const TotalScratch &);
  TotalScratch &operator=(const TotalScratch &);
};

static bool readTotalJobs(const string &path, vector<string> &jobs) {
  ifstream manifest(path.c_str());
  if (!manifest) {
    cerr << "Cannot open batch-total manifest: " << path << endl;
    return false;
  }
  string line;
  size_t number = 0;
  while (getline(manifest, line)) {
    ++number;
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty()) continue;
    if (line[0] != '/' || line.find_first_of("\t\r") != string::npos ||
        line.find('\0') != string::npos) {
      cerr << "Manifest line " << number
           << ": expected one absolute PDB path (no header or tabs)" << endl;
      return false;
    }
    struct stat info;
    if (stat(line.c_str(), &info) != 0 || !S_ISREG(info.st_mode)) {
      cerr << "Cannot read batch-total PDB: " << line << endl;
      return false;
    }
    ifstream input(line.c_str());
    if (!input) {
      cerr << "Cannot open batch-total PDB: " << line << endl;
      return false;
    }
    jobs.push_back(line);
  }
  if (!manifest.eof() || jobs.empty()) {
    cerr << "Batch-total manifest is empty or unreadable" << endl;
    return false;
  }
  return true;
}

static string csvQuoted(const string &value) {
  string result = "\"";
  for (size_t i = 0; i < value.size(); ++i) {
    if (value[i] == '"') result += '"';
    result += value[i];
  }
  return result + "\"";
}

static int runBatchTotal(int argc, char **argv) {
  string engine = "incremental";
  bool profiling = false;
  bool resume = false, verifyExisting = false;
  size_t flushEvery=16;
  if (argc < 4) {
    cerr << "Usage: BRiQ_Energy --batch-total JOBS.tsv TOTALS.csv [--engine legacy|shared|score|incremental] [--profile] [--resume|--resume-verify] [--flush-every N]" << endl;
    return 1;
  }
  for (int arg = 4; arg < argc; ++arg) {
    string option(argv[arg]);
    if (option == "--profile") profiling = true;
    else if (option == "--resume") resume = true;
    else if (option == "--resume-verify") {resume=true;verifyExisting=true;}
    else if (option == "--flush-every" && arg+1 < argc) {
      string value(argv[++arg]);
      if (value.empty() || value.find_first_not_of("0123456789")!=string::npos) {
        cerr<<"--flush-every requires an integer in 1..1024"<<endl;return 1;
      }
      char *end=NULL;errno=0;unsigned long parsed=strtoul(value.c_str(),&end,10);
      if (errno || *end || parsed<1 || parsed>1024) {
        cerr<<"--flush-every requires an integer in 1..1024"<<endl;return 1;
      }
      flushEvery=parsed;
    }
    else if (option == "--engine" && arg+1 < argc) engine = argv[++arg];
    else { cerr << "Unknown batch-total option: " << option << endl; return 1; }
  }
  if (engine != "legacy" && engine != "shared" && engine != "score" && engine != "incremental") {
    cerr << "Engine must be legacy, shared, score or incremental" << endl;
    return 1;
  }
  vector<string> jobs;
  if (!readTotalJobs(argv[2], jobs)) return 1;
  string output(argv[3]);
  if (output.empty() || output[0] != '/') {
    cerr << "Totals output must be an absolute container path" << endl;
    return 1;
  }
  const string checkpointPath=output+".state";
  const string identity=checkpointText(jobs,engine);
  vector<double> completed;
  struct stat outputInfo;
  const bool exists=lstat(output.c_str(),&outputInfo)==0;
  if (exists && (!resume || !S_ISREG(outputInfo.st_mode))) {
    cerr << "Totals CSV exists; use --resume for a verified regular checkpoint" << endl;
    return 1;
  }
  bool needsCheckpoint=true;
  size_t committedPrefix=0;
  vector<uint64_t> existingHashes;
  uint64_t csvHash=14695981039346656037ULL;
  struct stat checkpointInfo;
  if (!exists && lstat(checkpointPath.c_str(),&checkpointInfo)==0) {
    cerr << "Checkpoint exists without its totals CSV; choose a new output" << endl;
    return 1;
  }
  if (exists) {
    completed=readCompletedTotals(output,jobs);
    existingHashes=csvPrefixHashes(output);
    csvHash=existingHashes.back();
    ifstream state(checkpointPath.c_str());
    if (state) {
      ostringstream content;content<<state.rdbuf();
      committedPrefix=verifiedCheckpointPrefix(content.str(),identity,existingHashes);
      needsCheckpoint=false;
      verifyExisting=committedPrefix<completed.size();
    } else if (!verifyExisting) {
      cerr << "Legacy CSV has no checkpoint; use --resume-verify to recompute and verify its prefix" << endl;
      return 1;
    }
  }
  // Do not modify legacy data or create a sidecar until verification succeeds.
  int descriptor=open(output.c_str(),exists ? O_WRONLY|O_APPEND : O_WRONLY|O_CREAT|O_EXCL,0644);
  if (descriptor<0) {
    cerr << "Cannot create totals CSV (parent must exist; output must not exist): " << output << endl;
    return 1;
  }
  FILE *opened=fdopen(descriptor,exists ? "a" : "w");
  if (opened==NULL) {close(descriptor);return 1;}
  vector<char> csvBuffer(65536);
  unique_ptr<FILE,int (*)(FILE*)> totals(opened,fclose);
  if (setvbuf(totals.get(),csvBuffer.data(),_IOFBF,csvBuffer.size())!=0) return 1;
  if (!exists) {
    createCheckpoint(checkpointPath,identity);
    needsCheckpoint=false;
    const string header="job_index,input_pdb,E_total\n";
    if (fputs(header.c_str(),totals.get())<0 || fflush(totals.get())!=0) return 1;
    csvHash=csvHashUpdate(csvHash,header);
    recordCheckpoint(checkpointPath,0,csvHash);
  }
  if (exists && committedPrefix==jobs.size()) {
    cerr << "[resume] All jobs already verified by checkpoint" << endl;
    return 0;
  }
  cerr << "[resume] completed="<<completed.size()<<" verify-prefix="<<verifyExisting<<endl;
  cerr << "[batch-total] source_version=" << BRIQ_SOURCE_VERSION
       << " seed=0 weights=1,1,1" << endl;
  srand(0);
  cerr << "[batch-total] Loading shared energy tables" << endl;
  const auto resourceStart = chrono::steady_clock::now();
  EnergyTable tables;
  unique_ptr<BaseRotamerLib> sharedRotamers;
  unique_ptr<FragmentLibrary> sharedFragments;
  if (engine != "legacy") sharedRotamers.reset(new BaseRotamerLib());
  if (engine == "shared") sharedFragments.reset(new FragmentLibrary(sharedRotamers.get()));
  cerr << "[profile-resources] engine=" << engine << " seconds="
       << chrono::duration<double>(chrono::steady_clock::now()-resourceStart).count() << endl;
  IncrementalEnergy incremental;
  const size_t startIndex=committedPrefix;
  for (size_t index = startIndex; index < jobs.size(); ++index) {
    cerr << "[batch-total " << index + 1 << "/" << jobs.size() << "] "
         << jobs[index] << endl;
    const auto preparationStart = chrono::steady_clock::now();
    ScoringProfile profile;
    TotalScratch scratch;
    const bool light = engine=="score" || engine=="incremental";
    if (symlink(jobs[index].c_str(), scratch.pdb.c_str()) != 0 ||
        !(light ? prepareScoringInput(scratch.pdb,scratch.prepared)
                : prepareEnergyInput(scratch.pdb,scratch.prepared))) {
      cerr << "Cannot prepare batch-total job " << index << endl;
      return 1;
    }
    const double preparationSeconds = chrono::duration<double>(chrono::steady_clock::now()-preparationStart).count();
    vector<string> arguments = {
      "BRiQ_Energy", scratch.prepared,
      "--seed", "0", "--weights", "1", "1", "1"
    };
    vector<char *> pointers;
    for (size_t i = 0; i < arguments.size(); ++i) pointers.push_back(&arguments[i][0]);
    pointers.push_back(NULL);
    double energy = numeric_limits<double>::quiet_NaN();
    int status = runEnergy(static_cast<int>(arguments.size()), pointers.data(),
                           &tables, &energy, sharedRotamers.get(), sharedFragments.get(), &profile, light,
                           engine=="incremental" ? &incremental : NULL);
    if (profiling && engine=="incremental") cerr << "[incremental] " << index << ","
        << incremental.changedNodes << "," << incremental.pairEvaluations << ","
        << incremental.directedEvaluations << endl;
    if (profiling) cerr << "[profile] " << index << "," << preparationSeconds << ","
       << profile.rotamerSeconds << "," << profile.fragmentSeconds << ","
       << profile.treeSeconds << "," << profile.scoreSeconds << endl;
    if (status != 0 || !isfinite(energy)) {
      cerr << "Batch-total stopped at job " << index << endl;
      return status != 0 ? status : 1;
    }
    if (index<completed.size()) {
      if (fabs(energy-completed[index])>1e-9) {
        cerr << "Existing energy mismatch at job "<<index<<": "<<setprecision(17)
             <<completed[index]<<" vs "<<energy<<"; no new rows appended"<<endl;
        return 1;
      }
      if (index+1==completed.size() && needsCheckpoint) {
        createCheckpoint(checkpointPath,identity);
        needsCheckpoint=false;
      }
      if (index+1==completed.size()) recordCheckpoint(checkpointPath,completed.size(),csvHash);
      continue;
    }
    if (needsCheckpoint) {createCheckpoint(checkpointPath,identity);needsCheckpoint=false;}
    string quoted = csvQuoted(jobs[index]);
    ostringstream row;
    row << index << "," << quoted << "," << setprecision(17) << energy << "\n";
    string text=row.str();
    if (fwrite(text.data(),1,text.size(),totals.get())!=text.size()) {
      cerr << "Failed while writing totals CSV" << endl;
      return 1;
    }
    csvHash=csvHashUpdate(csvHash,text);
    if ((index+1)%flushEvery==0 || index+1==jobs.size()) {
      if (fflush(totals.get())!=0) {cerr<<"Cannot flush totals CSV"<<endl;return 1;}
      recordCheckpoint(checkpointPath,index+1,csvHash);
    }

  }
  if (fclose(totals.release()) != 0) {
    cerr << "Failed while closing totals CSV" << endl;
    return 1;
  }
  return 0;
}

#include "MutationBatch.h"

int main(int argc, char **argv) {
  try {
    if (argc > 1 && string(argv[1]) == "--mutate-total") return runMutationTotal(argc, argv);
    if (argc > 1 && string(argv[1]) == "--batch-total") return runBatchTotal(argc, argv);
    if (argc > 1 && string(argv[1]) == "--batch") return runBatch(argc, argv);
    return runEnergy(argc, argv);
  } catch (const exception &error) {
    cerr << "BRiQ_Energy failed: " << error.what() << endl;
    return 1;
  }
}
