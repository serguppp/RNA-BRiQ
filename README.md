# RNA-BRiQ

This repository contains codes for the following work: Peng Xiong, Ruibo Wu, Jian Zhan and Yaoqi Zhou, Robust RNA structure refinement by a nucleobase-centric sampling algorithm coupled with a backbone rotameric and quantum-mechanical-energy-scaled base-base knowledge-based potential. 

[![DOI](https://zenodo.org/badge/339727362.svg)](https://zenodo.org/badge/latestdoi/339727362)



This is free software for non-commercial users. If you have any questions, please contact Peng Xiong via email pengx_a@163.com


## Installation

### Prerequisites

Make sure these following prerequisites have been installed:
* C++ compiler: Needs to be C++11 compatible. Gcc 9.3.0 and Clang 9.1.0 tested.
* GNU Make
* Git
* CMake
* Wget and Tar - Optional: to download and extract the preprocessed data.

### Compilation

1. Clone from the repository:

```
git clone https://github.com/Jian-Zhan/RNA-BRiQ RNA-BRiQ
```

2. Change directory to `RNA-BRiQ/build/` and compile the codes:
```
mkdir RNA-BRiQ/build/
cd RNA-BRiQ/build/
cmake ../
make
```

### Data

Preprocessed data for the RNA-BRiQ can be downloaded from [http://servers.sparks-lab.org/downloads/RNA-BRiQ-data.tar.gz](http://servers.sparks-lab.org/downloads/RNA-BRiQ-data.tar.gz). The data should be downloaded, and extracted to `$BRiQ_DATAPATH` before continuing.

```
wget http://servers.sparks-lab.org/downloads/RNA-BRiQ-data.tar.gz

export BRiQ_DATAPATH=FILEPATH/BRiQ_data     ## Change "FILEPATH" to the real path to contain the preprocessed data
tar -xvz RNA-BRiQ-data.tar.gz --one-top-level=$BRiQ_DATAPATH
```

## Running the program

### Environment variables

The path containing BRiQ executables should be appended to `$PATH`. And `$BRiQ_DATAPATH` should be correctly set.

```
export BRiQ_DATAPATH=FILEPATH/BRiQ_data    ## Change "FILEPATH" to the real path containing the preprocessed data
export PATH=$PATH:FILEPATH/BRiQ/build/bin  ## Change "FILEPATH" to the real path containing the compiled codes
```

### Input files

Before running the program, you need to prepare the input file. An example files can be found at `demo/gcaa/input`

```
pdb gcaa.pdb    ## Initial structure in PDB format
seq GCGCAAGC    ## RNA nucleotide sequence
sec ((....))    ## Watson-Crick pairs
nwc ..(..)..    ## Non-Watson-Crick pairs
fixed 0 1 6 7   ## Index of the fixed nucleotides during structure sampling, starting from 0
```

If there is only RNA sequence but no reference PDB structure, `BRiQ_InitPDB` can generate an initial PDB structure from sequence. For example, the following command generates an initial PDB structure for RNA sequence `GCGCAAGC`.

```
BRiQ_InitPDB GCGCAAGC init.pdb
```

And the base pairing information can be extracted by `BRiQ_AssignSS` from the input PDB structure `$INPUTPDB` and stored to `$OUTFILE`:

```
BRiQ_AssignSS $INPUTPDB $OUTFILE
```

### RNA structure refinement

`BRiQ_Refinement` is the program to refine the initial structure of a RNA.

```
BRiQ_Refinement $INPUT $OUTPDB $RANDOMSEED
```

Where `$INPUT` is the input file, `$OUTPDB` is the output structure, and `$RANDOMSEED` is the random seed.

### RNA structure prediction

`BRiQ_Predict` is the program to predict the structure of a RNA.
```
BRiQ_Predict $INPUT $OUTPDB $RANDOMSEED
```

Where `$INPUT` is the input file, `$OUTPDB` is the output structure, and `$RANDOMSEED` is the random seed.

### Energy evaluation and pairwise components

Prepare the input format used by `BRiQ_Energy` directly from a PDB file:

```
BRiQ_Energy --prepare input.pdb input.briq
```

Calculate the global energy as before. To export one row per unordered residue
pair and a separate per-residue file for all remaining terms:

```
BRiQ_Energy input.briq --pair-energies pairs.csv
```

The command also writes `pairs.csv.nonpair.csv` unless a separate path is
given with `--nonpair`. The pair CSV columns are
`chain_i,residue_i,chain_j,residue_j,sep,contact_type,E_bb,E_bo,E_oo,E_clash,E_pair`.
`sep` is BRiQ's existing sequence-separation class (`0`, `1`, `2`, or `long`).
`contact_type` is determined independently of energy: WC/non-WC assignments
come from `AssignRNASS`; stacking requires adjacent connected residues with
base-normal angle <=30 degrees, plane distance <=4.0 A and base-center distance
<=8.0 A; base-ribose and base-phosphate labels use heavy-atom distances <=4.0 A.
`base_contact` uses a base-atom distance <=4.5 A. `unannotated` means none of
these explicit rules matched; it does not assert that the energy is zero or
that no physical contact exists.

The default `--pair-ownership pdb` assigns outgoing phosphate atoms to their
PDB owner (the next residue), matching BRiQ's own PDB writer. Earlier exports
grouped them by their storage node instead. Use `--pair-ownership node` to
reproduce that grouping; the CSV schema remains unchanged. In PDB mode `sep`
describes the reported endpoint residues, not necessarily the original node
separation used for each phosphate subterm.

`E_bb` is the six-dimensional base-base term; `E_bo` sums weighted base-ribose
and base-phosphate oxygen-table terms; `E_oo` sums ribose-ribose and
ribose-phosphate hydrogen-bond terms; `E_clash` sums the BRiQ atom-clash terms
from all pair interaction helpers. `E_pair` is the sum of these four columns.
The nonpair CSV has separate columns for same-residue base-ribose,
base-phosphate and ribose-phosphate terms, rotamer energy, and local phosphate
energy. Thus `sum(E_pair) + sum(E_nonpair)` reproduces the global score.
Local phosphate energy depends on an adjacent-ribose link and is anchored at
its upstream node for bookkeeping; it is not a physically one-body term.

Weights keep the existing defaults and can be set explicitly for comparison:

```
BRiQ_Energy input.briq --pair-energies pairs.csv --weights 0.5 1.25 0.75 --seed 0
```

The three values are the existing break, connect and clash weights. Original
BRiQ randomly selects loop closure flags; use `--seed` to make unequal
break/connect weighting repeatable. These flags are distinct from physical
input `break` positions. The default seed remains time-based.
The one-argument invocation and raw `--components COMPONENTS.tsv` output remain supported. Verify
the export, including weighted and two-chain cases, with:

```
python3 tests/test_pair_energy_export.py build/bin/BRiQ_Energy path/to/1EHZ.pdb
```

For a Docker-built image, run the same checks from the host with:

```
python3 tests/test_pair_energy_export.py --docker briq:tzok path/to/1EHZ.pdb
```

For independent validation of every local component against a separately
compiled, unmodified BRiQ reference, build an image from the current checkout
and run:

```
python tests/test_pair_energy_validation.py --docker briq:validated --build-reference
```

This uses the bundled 1RNA/1EHZ fixtures and absolute tolerance `1e-9`. It
retains logs, original directional helper outputs and all CSVs in the printed
artifact directory. On Windows, generated inputs use LF because the original
parser does not strip CRLF line endings. See the [validation report](docs/pair-energy-validation.md)
for term definitions, ownership conventions, results and limitations.

### Demonstrations

* [demo/refinement/](demo/refinement/): RNA structure refinement.
* [demo/partial_predict/](demo/partial_predict/): RNA structure prediction from an initial reference structure.
* [demo/predict_from_seq_SS/](demo/predict_from_seq_SS/): RNA structure prediction from sequence and pre-assigned secondary structure information.
