# RNA-BRiQ — energy exporter

Export RNA structure energies with BRiQ, including per-pair components and the
remaining backbone and local terms.

This fork builds on [RNA-BRiQ](https://github.com/Jian-Zhan/RNA-BRiQ) and
[Tomasz Żok's energy scorer](https://github.com/tzok/RNA-BRiQ).
It preserves the original scoring formulas and adds CSV export and validation.
It does not calculate frustration indices.

## Quick start

With `BRiQ_Energy` installed and `BRiQ_DATAPATH` configured:

```bash
BRiQ_Energy molecule.pdb --output-dir results --seed 0
```

The output directory must be new or empty; its parent must already exist.
A successful run produces exactly these five files:

| File | Contents |
| --- | --- |
| `molecule.pdb` | Unchanged copy of the input structure |
| `molecule.briq` | Generated sequence, pairing annotations and chain breaks |
| `pairs.csv` | Energy components for every unordered residue pair |
| `nonpair.csv` | Same-residue terms, rotamer energies and local phosphate energies |
| `run.json` | Program and source versions, scoring options, random seed and total energy |

The PDB keeps its original filename; the `.briq` file uses the same stem.
For example, `1EHZ.pdb` produces `1EHZ.pdb` and `1EHZ.briq`.
The generated configuration refers to the copied PDB by a relative path, so the
directory can be moved. Run later scoring commands from inside that directory.
Diagnostics go to the console. A failed run may leave partial output.

### Docker on Windows

If you have the locally built `briq:output-dir` image, run this in PowerShell.
Change the input path for each structure:

```powershell
$pdb = (Resolve-Path '.\molecule.pdb').Path
$name = [System.IO.Path]::GetFileName($pdb)
$stem = [System.IO.Path]::GetFileNameWithoutExtension($pdb)
$out = [System.IO.Path]::GetFullPath("$PWD\results\$stem")

New-Item -ItemType Directory -Force -Path $out | Out-Null

docker run --rm `
  --mount "type=bind,source=$pdb,target=/data/$name,readonly" `
  --mount "type=bind,source=$out,target=/output" `
  briq:output-dir `
  BRiQ_Energy "/data/$name" --output-dir /output --seed 0
```

`briq:output-dir` is a local image tag, not a published container image.

## Build from source

Use Linux or WSL with Git, CMake, GNU Make and a C++11 compiler.
The energy exporter also requires the BRiQ parameter tables.

```bash
git clone https://github.com/serguppp/RNA-BRiQ.git
cd RNA-BRiQ
mkdir build
cd build
cmake ..
make -j2 BRiQ_Energy
cd ..
export PATH="$PWD/build/bin:$PATH"
```

Download and extract the [upstream parameter archive](http://servers.sparks-lab.org/downloads/RNA-BRiQ-data.tar.gz).
Set `BRiQ_DATAPATH` to the extracted parameter directory:

```bash
export BRiQ_DATAPATH=/absolute/path/to/BRiQ_data
```

For reproducible results, keep the executable and parameter tables together,
or preserve the Docker image. `run.json` records scoring options, but does not
embed or checksum the parameter tables.

## Scoring options

```bash
BRiQ_Energy molecule.pdb --output-dir results \
  --weights 1 1 1 --pair-ownership pdb --seed 0
```

| Option | Meaning | Default |
| --- | --- | --- |
| `--weights BREAK CONNECT CLASH` | Three outer BRiQ scoring weights | `1 1 1` |
| `--pair-ownership pdb\|node` | Assign phosphates to their PDB residue or internal storage node | `pdb` |
| `--seed N` | Fix the random seed for repeatable runs | Time-based |

Use a fixed seed when comparing structures or scoring options. BRiQ uses
random loop-closure flags that can affect scores when BREAK and CONNECT differ.

### Reading the energies

`pairs.csv` contains **all residue pairs**, including pairs with zero energy.

| Column | Energy contribution |
| --- | --- |
| `E_bb` | Six-dimensional base–base potential |
| `E_bo` | Base–ribose and base–phosphate oxygen-table terms |
| `E_oo` | Ribose–ribose and ribose–phosphate hydrogen-bond terms |
| `E_clash` | Pairwise steric terms |
| `E_pair` | Sum of the four components above |

The complete score is recovered, within floating-point precision, by:

```text
total energy = sum(pairs.csv:E_pair) + sum(nonpair.csv:E_nonpair)
```

`nonpair.csv` groups remaining terms by residue for bookkeeping. Local
phosphate energy depends on a link between adjacent riboses, so it is not a
single-residue physical interaction.

`contact_type` provides geometric and BRiQ pairing annotations. These labels
are not a complete contact classification; `unannotated` does not mean there
is no contact, and zero energy does not mean there is no contact.

See the [technical report](docs/pair-energy-validation.md) for term definitions,
phosphate ownership, sequence-separation conventions and limitations.

## Advanced export

To prepare and score separately, or choose individual output paths:

```bash
BRiQ_Energy --prepare molecule.pdb input.briq
BRiQ_Energy input.briq --seed 0 \
  --pair-energies pairs.csv --nonpair nonpair.csv
```

Score without CSV export:

```bash
BRiQ_Energy input.briq --seed 0
```

Add `--components components.tsv` for raw directional helper terms.
Individual output flags cannot be combined with `--output-dir`.
If `--nonpair` is omitted during pair export, its filename defaults to
`<pair-output-path>.nonpair.csv`.

## Validation

Run from the repository root with Docker and Python installed:

```bash
# Five-file output, score reconstruction and overwrite protection
python tests/test_output_directory.py --docker briq:output-dir

# Compare exported components against a separately built original BRiQ scorer
python tests/test_pair_energy_validation.py --docker briq:output-dir --build-reference
```

The tests use the bundled 1RNA and 1EHZ structures. Detailed results and known
limitations are in the [validation report](docs/pair-energy-validation.md).
Agreement with original BRiQ validates the export; it does not establish the
scientific accuracy of the potential or a frustration calculation.

## Other BRiQ tools

Build all targets with `make -j2` from `build/` to use the original prediction
and refinement tools. Examples are in
[demo/refinement/](demo/refinement/),
[demo/partial_predict/](demo/partial_predict/) and
[demo/predict_from_seq_SS/](demo/predict_from_seq_SS/).

## Reference and usage terms

Peng Xiong, Ruibo Wu, Jian Zhan and Yaoqi Zhou.
*Robust RNA structure refinement by a nucleobase-centric sampling algorithm
coupled with a backbone rotameric and quantum-mechanical-energy-scaled
base-base knowledge-based potential.*
[Nature Communications (2021)](https://doi.org/10.1038/s41467-021-23100-4).

Upstream states that the software is free for non-commercial users.
For questions about upstream usage terms, contact Peng Xiong at
pengx_a@163.com.
