# Scoring performance and correctness

The total-only batch scorer reuses potential and rotamer tables, skips sampling
setup, and optionally caches unchanged energy contributions. It still computes
the complete BRiQ total: mutations affect interactions with the environment,
not just the energy of the selected residue pair.

## Measurements on this workstation

Windows host, Docker Desktop Linux engine, input PDBs on a Windows bind mount.
Only one scoring process ran at a time. Both engines below include table loading,
input fingerprinting, preparation, checkpoint/CSV output and container overhead;
mutants were prepared beforehand. Full scoring ran before incremental scoring,
so filesystem cache/load differences also affect wall-time ratios.

| Structure | N | Contacts C | Variants | Full scorer | Incremental scorer |
| --- | ---: | ---: | ---: | ---: | ---: |
| 1RNA | 28 | 60 | 960 | 85.36 s | 64.83 s |
| 1EHZ | 76 | 173 | 2768 | 256.46 s | 174.25 s |

Validation details: [integration-test output](performance-validation.txt).
The corrected native full-export totals were -28.606109849613684 for 1RNA and
-136.42678705881318 for 1EHZ, with seed 0 and outer weights 1/1/1. These values
refer to the bundled fixtures and this parameter-table image.

The maximum full-vs-incremental energy difference was zero for all 3728 variants
in separate container runs of the corrected build. This validates agreement of
these implementations, not biological validity of a frustration index. The
observed batch wall-time ratios are 1.32x and 1.47x; these do not include the much
larger gain from eliminating per-variant fragment-library loading earlier.

On a mixed 16-job sample the corrected legacy engine took 90.14 s, shared
libraries 37.51 s, scoring-only 31.93 s and incremental 33.96 s. Initialization
dominates this small sample; incremental caching is not guaranteed to help
every workload. The largest earlier cost was fragment-library loading on every
job (about 3.1-3.6 s per job in that sample). Do not extrapolate a universal
speedup factor from these runs. See [engine comparison](performance-engines.json).

Raw measurements: [final scoring](performance-final.json),
[sequential Python generation](performance-generation.json).
Table initialization alone took about 30-32 s per batch. The unchanged serial
Python generator took 29.23 s for 1RNA and 237.95 s for 1EHZ in an earlier timed
run on the same workstation. Generation plus incremental scoring therefore sums
to about 94 s and 412 s. Separate native-export regression runs took 36.27 s
for 1RNA and 37.12 s for 1EHZ. Adding these gives approximately **2-2.5 minutes
for 1RNA** and **7-8 minutes for
1EHZ**, respectively. These are assembled end-to-end estimates, not timed single
CLI invocations. No parallel Python generation is required for those estimates.

The energy-evaluation portion of the incremental batch profiles sums to only
about 0.51 s (1RNA) and 2.68 s (1EHZ). Parsing, geometry/phosphate construction,
initialization, input hashing and filesystem/output overhead account for the
rest; optimizing the potential evaluator alone cannot remove that work.

## Correctness issue discovered during benchmarking

The original phosphate builder indexed a vector containing one error value per
level-1 library cluster with a dihedral-grid index. AddressSanitizer reproduced
an invalid read (reported as heap-use-after-free). Allocator-dependent values
could select a different phosphate refinement search and produce different
energies across otherwise equivalent runs. The corrected code indexes by the
selected level-1 cluster and uses bounds checking. The standalone sanitizer
probe exercises both improper-angle signs on 1RNA and 1EHZ (54 and 150 probes).

**Regenerate the native score and every mutant with the corrected version.**
Do not mix scores from before and after this fix. Historical stage benchmark
JSONs record execution times, not a trustworthy scientific reference for the
old energies. Earlier same-build agreement did not reveal this invalid read.

## Implemented changes

- Potential and rotamer libraries are loaded once per total-only batch.
- Scoring-only construction skips simulation connections, fragment libraries,
  sampling caches and duplicate energy evaluations.
- Total-only preparation omits secondary-structure annotation, retaining the
  same residue sequence and physical/chain connectivity.
- Incremental scoring caches contributions only for exactly unchanged node
  coordinates/types and separation; a changed node invalidates all its
  interactions. Phosphate geometry changes are included in this comparison.
- CSV output is buffered, with a checkpoint every 16 rows by default.
- Resume checks source version, scoring engine, parameters, manifest/file
  contents and committed CSV checksums before skipping a completed prefix.
  An uncommitted tail is recomputed and compared before appending.
- Successful total-only jobs suppress per-node diagnostics. Errors retain them.

Default engine: `incremental`. Use `--engine score` for full recomputation,
`--engine shared` to retain sampling construction with shared libraries, or
`--engine legacy` to reproduce per-job resource loading for profiling.

## Historical validation scripts

`tests/test_pair_energy_validation.py` retains a pinned upstream oracle and
frozen pre-fix values for the historical export study. Its unmodified upstream
phosphate builder contains the invalid read, so it is not an oracle for the
corrected model preparation. Current regression checks compare full export,
full scoring and incremental scoring within one corrected version, including
independent processes and repeated inputs. Five-file reconstruction checks
remain active; historical hardcoded totals were removed from those mode tests.

## Complexity

Let N be the number of retained nucleotides and C the number of selected
contacts. There are V = 16C requested variants, including the native sequence
combination according to the chosen geometry preparation protocol.

The current Python generator repeatedly parses and copies entire target chains
and writes full PDBs. Copying and file I/O alone require O(CN) work, with
additional ModeRNA preparation, sorting and format-conversion overhead.
Mutation-specific work can vary.
Full scoring visits O(N^2) residue pairs per variant, giving O(CN^2) work.
Incremental scoring evaluates roughly O(kN) expensive interactions when k nodes
change, but still traverses O(N^2) cached entries. It therefore retains the same
worst-case asymptotic complexity. Per-variant construction also rebuilds the
quadratic sequence-separation table, and the cache copies that table. These
steps must be eliminated/reused as well to achieve linear per-variant work.
Rebuilding coordinates and phosphate groups
also occurs for every variant. Dense changes invalidate the whole cache.

If C grows proportionally to N, full enumeration is approximately cubic.
If all possible residue pairs are selected, C = O(N^2) and the worst case is
quartic. Neither relationship determines wall time for a particular structure.
The native all-pairs CSV and scorer caches require O(N^2) memory; retained full
mutant PDBs require O(CN) disk space. Shared potential tables add a large fixed
memory requirement independent of N.

## Estimating another structure

Measure a spread of K contacts in the actual structure, with all 16 variants per
contact, the same geometry protocol and one persistent batch. Separate fixed
startup/native-scoring cost T0 from the variable generation/scoring cost:

    T_estimated = T0 + (C / K) * (T_pilot - T0)

Use several contacts from different regions and report a range. Do not scale a
short pilot by C/K without subtracting initialization: loading the tables once
can dominate the pilot. Profiling reports resource initialization separately.
A more general empirical model is:

    T = T0 + 16 C (a N + b N^2 + d)

Fit its constants on the same machine, image, filesystem and scoring engine.
Two structures alone do not provide a reliable fitted scaling law. An N^3
estimate assumes similar contact density and a dominating quadratic scorer;
those assumptions do not describe every stage of the optimized pipeline.

## Further sequential optimizations

1. Parse and prepare the immutable ModeRNA/Gemmi structure once. Precompute four
   replacement residues per position and combine them for each contact, instead
   of repeatedly deep-copying and rebuilding entire chains. Validate coordinates,
   residue identities, occupancies and scores against the existing protocol.
2. Replace full PDB exchange with a persistent in-memory scoring interface:
   initialize the native structure once, send mutated residue data/coordinates,
   rebuild affected phosphates, evaluate, and restore the baseline.
3. Update cached interactions using an explicit list of changed nodes. Avoid
   scanning all N^2 entries merely to discover they are unchanged. Reuse the
   separation/connectivity table when the backbone topology is unchanged. Maintain a
   baseline total plus affected-term differences; allow and test small changes
   from floating-point summation order.
4. Precompute the four-state node data and reusable interaction tables where
   geometry independence has been demonstrated. Include adjacent phosphates
   and all environment terms; do not reduce the calculation to E_pair(i,j).
5. Add spatial neighbor lists only after auditing the exact support/cutoffs of
   every energy term. Dropping a nonzero term would change the energy model.
6. Obtain native annotation/export and mutant totals in one long-lived process
   to avoid loading the large tables in two separate containers.
7. Reduce Windows bind-mount I/O by staging inputs inside the Linux filesystem
   or transmitting a compact archive. Preserve provenance and benchmark the
   staging cost as part of wall time.
8. Deduplicate only genuinely identical prepared geometries plus all scoring
   parameters. With a demonstrated, context-independent four-state geometry per
   site and one common reconstructed baseline, 16C requests can collapse to at
   most 1 + 3N_active + 9C distinct sequence states (one baseline, single changes,
   and double changes). This bound is not automatically valid for the current
   contact-specific reconstruction. A native-sequence substitution is not
   automatically identical to an untouched experimental PDB.

These are implementation ideas, not measured additional speedups. Checkpointing
saves repeated work after interruption; it is not an asymptotic scoring speedup.
Changing to a coarse-grained potential is a separate model comparison and needs
its own validation and benchmark.

## Reproduce tests

Run from the repository root with the built image and Python:

```powershell
python tests/test_resume.py --image briq:batch
python tests/benchmark_scoring.py --workspace C:/path/to/prepared/benchmark `
  --manifest jobs-benchmark.tsv --engines legacy shared score incremental
```

The benchmark asserts finite, ordered results and energy agreement within 1e-9.
Use fresh output names; the scorer refuses accidental overwrite.

For the lightweight phosphate memory-safety regression, inside the build image:

```sh
g++ -std=c++11 -g -O1 -fsanitize=address -fno-omit-frame-pointer -I. \
  tests/check_phosphate.cpp forcefield/src/PO3Builder.cpp \
  build-batch/lib/libbriq.a -pthread -o /tmp/check-phosphate
ASAN_OPTIONS=detect_leaks=0 /tmp/check-phosphate tests/fixtures/1RNA.pdb
ASAN_OPTIONS=detect_leaks=0 /tmp/check-phosphate tests/fixtures/1EHZ.pdb
```

The phosphate builder is instrumented directly; the existing static library
supplies its geometry/IO dependencies. This targeted test does not claim the
entire upstream application is free of memory defects.
