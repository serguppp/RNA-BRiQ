/*
 * RunEnergy.cpp
 *
 *  Created on: 2025-04-28
 *      Author: tzok
 */

#include "pred/BRFoldingTree.h"
#include "model/AssignRNASS.h"
#include "model/ProteinRep.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <sstream>
#include <set>
#include <vector>

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

void checkInputFile(const string &inputFile) {
  ifstream file;
  cout << "check input file: " << endl;
  file.open(inputFile.c_str(), ios::in);
  if (!file.is_open()) {
    cout << "fail to open inputFile: " << inputFile << endl;
    exit(0);
  }
  file.close();

  NSPtools::InputParser input(inputFile);
  string pdbFile = input.getValue("pdb");
  string baseSeq = input.getValue("seq");
  string baseSec = input.getValue("sec");
  string nwcSec = input.getValue("nwc");

  if (pdbFile == "") {
    cout << "invalid inputFile: initial pdb required" << endl;
    exit(0);
  }

  if (baseSeq == "") {
    cout << "invalid inputFile: base sequence required" << endl;
    exit(0);
  }

  if (baseSec == "") {
    cout << "invalid inputFile: Watson-Crick pairing sequence required" << endl;
    exit(0);
  }

  if (nwcSec == "") {
    cout << "invalid inputFile: non-Watson-Crick pairing sequence required"
         << endl;
    exit(0);
  }

  if (baseSec.length() != baseSeq.length()) {
    cout << "the length of Watson-Crick pairing sequence must equal to the "
            "length of base sequence"
         << endl;
    exit(0);
  }

  if (nwcSec.length() != baseSeq.length()) {
    cout << "the length of non-Watson-Crick pairing sequence must equal to the "
            "length of base sequence"
         << endl;
    exit(0);
  }

  file.open(pdbFile.c_str(), ios::in);
  if (!file.is_open()) {
    cout << "fail to open pdb file: " << pdbFile << endl;
    exit(0);
  }
}

void printHelp() {
  cout << "Usage: BRiQ_Energy INPUT_FILE [--pair-energies PAIRS.csv]"
          " [--nonpair NONPAIR.csv] [--components COMPONENTS.tsv]"
          " [--weights BREAK CONNECT CLASH] [--pair-ownership pdb|node]"
          " [--seed SEED]" << endl;
  cout << "       BRiQ_Energy --prepare INPUT.pdb OUTPUT_FILE" << endl;
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

int main(int argc, char **argv) {
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
    if (option == "--pair-energies") {
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

  checkInputFile(inputFile);

  EnergyTable *et = new EnergyTable();
  BRFoldingTree *ft = new BRFoldingTree(inputFile, et);
  double energy;
  ofstream components;
  if (!componentFile.empty()) {
    components.open(componentFile.c_str(), ios::out);
    if (!components.is_open()) {
      cerr << "fail to open component output file: " << componentFile << endl;
      delete et;
      delete ft;
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
      delete et;
      delete ft;
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
      delete et;
      delete ft;
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
        delete et;
        delete ft;
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
      delete et;
      delete ft;
      return 1;
    }
  } else if (components.is_open()) {
    energy = ft->totalEnergy(breakWeight, connectWeight, clashWeight, 0.0,
                             false, &components);
  } else {
    energy = ft->totalEnergy(breakWeight, connectWeight, clashWeight, 0.0,
                             false);
  }
  if (components.is_open()) {
    components.close();
    if (!components.good()) {
      cerr << "failed while writing component output file: " << componentFile
           << endl;
      delete et;
      delete ft;
      return 1;
    }
  }
  if (!isfinite(energy)) {
    cerr << "BRiQ produced a non-finite energy" << endl;
    delete ft;
    delete et;
    return 1;
  }
  if (exportPairs)
    cerr << "pair_ownership " << (pdbOwnership ? "pdb" : "node") << endl;
  cerr << inputFile << " " << setprecision(17) << energy << endl;

  delete et;
  delete ft;
}
