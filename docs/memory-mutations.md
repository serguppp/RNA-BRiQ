# Mutations in one scoring process

`--mutate-total` evaluates the reference and every A/C/G/U substitution of
each selected reference contact in one process. It loads EnergyTable and
BaseRotamerLib once and retains the native scoring tree. It does not generate
or parse a complete mutant PDB for each substitution.

```bash
BRiQ_Energy --mutate-total reference.pdb geometry.tsv output
```

The parent of `output` must exist. The output directory must be empty or absent.
The reference must be a single-model PDB with valid RNA residues. Blank chain
identifiers and PDB insertion codes are supported; duplicate chain/residue
identities are rejected. The geometry cache must not change reference chain
connectivity. This is an equal-weight, seed-zero scorer, matching the existing
total-batch defaults, with PDB phosphate ownership for native CSV exports.

For the separate Python application, see the [code to copy into briq.py and pipeline.py](python-frustratometer-integration.md).

## Geometry supplied by the caller

The TSV is one file with this exact header:

```text
chain<TAB>residue<TAB>base<TAB>atom<TAB>x<TAB>y<TAB>z
```

Every row contains one heavy atom of a cached canonical A/C/G/U residue, in
the global coordinate system of the reference. A cache entry is identified by
`(chain, residue, base)`; atom names within an entry must be unique. Supply all
four entries for every residue participating in selected contacts. Supplying
all reference residues allows contact discovery inside C++, with no preliminary
native-scoring invocation. An empty chain field denotes a blank PDB chain.

Generate the coordinates using the same mutation procedure as the PDB pipeline.
To reproduce PDB scoring, export three decimal places. The reader intentionally
uses the original PDB reader's float conversion before BRiQ standardization.
It does not change the template fitting, minimize geometry, or skip same-base
reconstructions. Base frames are computed after standardization, as in the file
scorer. Modified reference residues retain the existing upstream interpretation;
mutation entries themselves must be canonical. This mode does not resolve the
scientific choice of reference-versus-reconstructed geometry.

## Execution

1. Validate the reference and geometry cache, including connectivity.
2. Discover contacts with the existing reference annotation/classification.
3. Load scoring libraries once and export the reference with the retained tree.
4. Standardize each cached residue once and construct its BRNode/rotamer once.
5. For each contact, select two cached nodes for each of 16 combinations in
   `AA, AC, AG, AU, CA, ... UU` order. All other nodes remain reference nodes.
6. Update dependent phosphate links and evaluate the full total energy using
   the existing exact incremental scorer. Restore original node pointers.

Phosphate construction is also memoized for each ordered pair of immutable
reference/cached node geometries. Its inputs are the first node's cs2, the
second node's ribose coordinates, and their improper angles. Reusing this
function result changes neither geometry nor energy. Changed neighboring
phosphates are included in incremental invalidation. No interaction cutoff or
pair-only approximation is introduced.

## Outputs

The native five-file export is preserved:

```text
reference.pdb           # uses the actual input basename
reference.briq
pairs.csv
nonpair.csv
run.json
contact_variants.csv    # all 16 energies per selected reference contact
energy_components.csv  # molecular component sums: reference + all variants
mutation_run.json      # completion metadata and timings
```

`contact_variants.csv` columns are:

```text
chain_i,residue_i,chain_j,residue_j,contact_type,base_i,base_j,native_E_total,E_total,delta_E_total
```

`E_total` is the full molecular score. Contact classification and ordering match
native `pairs.csv` after removing `unannotated` rows. No frustration thresholds
or index definition are changed. Output is written as `.tmp` and finalized only
after scoring succeeds; `mutation_run.json` records a successful complete run.
This new mode deliberately refuses existing results and does not implement
resume. The existing `--batch-total` checkpoint/resume interface is unchanged.

The current mode also writes [molecule-wide component sums](energy-components.md)
without changing the contact-variant CSV schema or Python invocation.

## Validation

```bash
python tests/test_mutation_total.py --image briq:memory --verify-full
```

Bundled fixtures contain 1RNA and 1EHZ reference/cache inputs and independent
PDB-batch energies produced with corrected source d0a1f2b. Their SHA-256 hashes
and generation protocol are recorded in `tests/fixtures/memory/provenance.json`.
The test checks all 3,728 energies, all 16 combinations per contact, annotations,
reference energy and energy differences. It also rejects bad headers, empty
caches, repeated atoms, non-finite coordinates, unknown residues, missing variants,
changed connectivity and existing results. `--verify-full` independently rebuilds
phosphates and evaluates every energy term for each variant, bypassing both
phosphate memoization and the incremental energy cache for the comparison.

## Measured serial performance

Windows/Docker Desktop measurements on 2026-11-04, including Python cache
preparation/export and Docker execution, excluding interpreter startup:

| Structure | Contacts | Variants | Complete run | Library loading | Variant loop |
|---|---:|---:|---:|---:|---:|
| 1RNA | 60 | 960 | 29.60 s | 26.84 s | 0.16 s |
| 1EHZ | 173 | 2,768 | 31.91 s | 27.77 s | 0.80 s |

The previous measured 1EHZ cache/PDB pipeline took 251.12 s: approximately
7.9 times slower than this run. The 1RNA historical 2–2.5 minute estimate was
from an earlier pipeline and is not a directly matched before/after benchmark.
All energies and native pair/nonpair CSV rows match the current PDB pipeline
exactly for both structures. These are single-run wall times, not a timing
guarantee; parameter-table and filesystem caches affect initialization.

Memory holds cached nodes/rotamers and incremental pair tables; the latter remain
O(N²). The incremental scorer still scans the pair table for each variant, so
this does not claim a new asymptotic bound. Input/output volume falls from
O(CN) full-structure PDB data to O(N) cached geometry plus O(C) energy rows.
