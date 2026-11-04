# RNA-BRiQ — energy exporter

Export RNA structure energies with BRiQ, including per-pair components and the
remaining backbone and local terms.

This fork builds on [RNA-BRiQ](https://github.com/Jian-Zhan/RNA-BRiQ) and
[Tomasz Żok's energy scorer](https://github.com/tzok/RNA-BRiQ).
It retains the scoring formulas and adds CSV export and validation.
Phosphate construction includes a [memory-safety correction](docs/performance.md);
recompute native and mutant energies together when upgrading older runs.
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

After building the current `briq:batch` image (instructions below), run this in PowerShell.
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
  briq:batch `
  BRiQ_Energy "/data/$name" --output-dir /output --seed 0
```

`briq:batch` is a local image tag, not a published container image.

## Batch scoring

`BRiQ_Energy --batch JOBS.tsv` scores multiple PDB files in one process and
loads the energy tables once. Each structure still gets its own five-file output.
Batch mode uses seed **0**, weights **1/1/1**, and **PDB** phosphate ownership.
Single-structure scoring remains available with the existing options.

The manifest has no header. Each row contains an absolute **container** PDB path,
a tab character, and an absolute container output-directory path. Spaces are
allowed; tabs and newlines within paths are not. LF and CRLF files are accepted.
For example (the gap between the paths must be a tab):

```text
/work/structures/1RNA.pdb	/work/scores/original
/work/structures/1RNA_A1A.pdb	/work/scores/A1A
```

All inputs and output parents must exist. Output directories must be new or empty
and unique. The whole manifest is checked before tables are loaded or job outputs
are created. Jobs run sequentially, resetting the seed for each structure; scoring
stops on an error. A failed run may leave completed or partial job outputs.

### Build the batch image

From this repository in PowerShell, using the existing local image containing
BRiQ's parameter tables:

```powershell
$sourceVersion = git describe --always --dirty
docker build -f Dockerfile.batch `
  --build-arg "BRIQ_SOURCE_VERSION=$sourceVersion" `
  -t briq:batch .
```

The build reuses `briq:output-dir`. To use a different compatible local base image,
add `--build-arg BRIQ_BASE_IMAGE=YOUR_IMAGE`. The parameter tables are inherited;
they are not downloaded again. The source version is required so `run.json`
identifies the scorer build.

### Run a batch

Prepare `jobs.tsv` and the directories inside a workspace, then mount that whole
workspace. For example, if the manifest paths begin with `/work/`:

```powershell
$workspace = (Resolve-Path '.\batch-work').Path
docker run --rm `
  --mount "type=bind,source=$workspace,target=/work" `
  briq:batch BRiQ_Energy --batch /work/jobs.tsv
```

Compare separate processes with a batch using the same corrected image:

```powershell
python tests/test_batch_energy.py --docker briq:batch
python tests/test_output_directory.py --docker briq:batch
```

The regression compares total, pair and nonpair components for 1RNA and 1EHZ,
repeats an input after a different structure, and checks manifest preflight and
overwrite protection. It reports elapsed times for three separate runs versus
one three-job batch. This is a small benchmark, not a prediction of full-scan time.
Measured results are recorded in [batch-validation.md](docs/batch-validation.md).

## Total-only batch scoring

Use this mode when only the total energy of each variant is needed:

```text
BRiQ_Energy --batch-total JOBS.tsv TOTALS.csv
```

Unlike `--batch`, this manifest has **one absolute container PDB path per row**,
without a header or an output-directory column. Blank rows and CRLF are accepted.
Repeated inputs are allowed and produce separate rows in manifest order.

```text
/work/structures/variant_AA.pdb
/work/structures/variant_AC.pdb
```

```powershell
$workspace = (Resolve-Path '.\batch-work').Path
docker run --rm `
  --mount "type=bind,source=$workspace,target=/work" `
  briq:batch BRiQ_Energy --batch-total /work/jobs-total.tsv /work/totals.csv
```

`TOTALS.csv` must be an absolute container path with an existing parent; the file
must not already exist. All manifest paths are checked before creating the CSV
or loading potential tables. Output columns are `job_index,input_pdb,E_total`;
indices start at zero and paths are CSV-escaped. Connect contact IDs and mutation
identities to these indices in the caller instead of parsing PDB filenames.

One set of potential tables is loaded. Each job uses seed 0 and scoring weights
1/1/1, with equivalent sequence/connectivity preparation and total-energy terms.
Preparation uses a private temporary directory inside the container, removed
after each job. The mode does not create per-variant output directories, copy
PDBs into outputs, annotate exported contacts, or write pair/nonpair/run files.
It still prepares and scores each entire structure. Source version and fixed
parameters are written to stderr; retain that log for provenance.

The CSV is flushed every 16 jobs by default (`--flush-every 1` restores per-job
flushing). A `.state` sidecar records input fingerprints and checksums of flushed
CSV prefixes. On a runtime failure the CSV can contain only a prefix of results:
check the process exit status and expected row count. Existing results are never
overwritten; explicit resume options are described below. Existing full-export modes are retained
for native contact discovery and detailed inspection of selected variants.

```powershell
python tests/test_batch_total.py --docker briq:batch
```

### Engines and profiling

`--batch-total` accepts `--engine`:

| Engine | Behavior |
| --- | --- |
| `legacy` | Reload libraries and build the original sampling tree per job |
| `shared` | Share rotamer/fragment libraries while retaining the original tree |
| `score` | Skip sampling setup and evaluate every total-energy term |
| `incremental` (default) | Reuse terms whose geometry, identities and connectivity are unchanged |

All engines use seed 0 and outer weights 1/1/1. The incremental engine invalidates
all interactions of changed nodes, including their environment and neighboring
phosphates. It preserves total energy and accumulation order; it does not replace
the total score with the energy of only the mutated pair. Detailed native contact
discovery continues to use the full exporter.

Add `--profile` to log per-job preparation, library loading, tree setup and energy
calculation times. Resource initialization is logged separately. See the measured
[performance and validation report](docs/performance.md). A reusable benchmark is
`tests/benchmark_scoring.py --workspace WORKSPACE --image briq:batch`.

### Resume a total-only batch

```text
BRiQ_Energy --batch-total JOBS.tsv TOTALS.csv --resume
```

Keep the same image, parameter tables, input files, manifest and `.state` sidecar.
The sidecar checks source version, engine, fixed parameters, input contents and
committed CSV rows, using non-cryptographic FNV-1a checksums for accidental changes.
The default flush/checkpoint interval is 16 jobs (`--flush-every 1` writes each job).
Any valid uncommitted tail is recomputed and compared before new rows are appended.
A completed verified batch exits without loading tables. Invalid/truncated CSVs,
changed inputs or altered committed energies are rejected without silent repair.

For older CSVs without `.state`, use `--resume-verify`. It recomputes and checks
**all existing energies** within `1e-9` before creating a checkpoint or appending.
The existing prefix is preserved; incompatible scores stop the run.

```powershell
python tests/test_resume.py --image briq:batch
```

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

## Mutational analysis without mutant PDB files

Use `BRiQ_Energy --mutate-total reference.pdb geometry.tsv output` to score
native contacts and all their 16 A/C/G/U substitutions in one process.
The caller supplies one cache of residue geometries; BRiQ reuses its potential,
node and phosphate calculations. It returns full molecular energies and the
native export, plus [component sums for the reference and every variant](docs/energy-components.md).
See [input format, results and measured timings](docs/memory-mutations.md).

## Validation

The corrected scorer fixes an invalid phosphate-library read that could make
older scores depend on memory layout. Recompute native and mutant scores with
one corrected image; do not mix them with historical results. See the
[correctness and performance report](docs/performance.md).


Run from the repository root with Docker and Python installed:

```bash
# Five-file output, score reconstruction and overwrite protection
python tests/test_output_directory.py --docker briq:batch

# Current full-export/total-only agreement and checkpoint tests
python tests/test_batch_total.py --docker briq:batch
python tests/test_resume.py --image briq:batch
```

`test_pair_energy_validation.py` and its frozen reference values document the
historical exporter validation against unmodified upstream BRiQ. That upstream
reference contains the phosphate bug and is not an oracle for corrected scores.

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
