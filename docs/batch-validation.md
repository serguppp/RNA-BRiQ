# Batch scoring validation

Validated on 2026-11-03 using Docker Desktop's Linux engine.

## Numerical compatibility

`tests/test_batch_energy.py --docker briq:batch --reference briq:output-dir`
passed for the sequence **1RNA, 1EHZ, 1RNA**. The previous image was used as
an independent executable reference. All total energies and every numeric
component in `pairs.csv` and `nonpair.csv` matched with absolute tolerance
`1e-9`. Pair identities, annotations and copied inputs were also checked.
Repeating 1RNA after 1EHZ produced byte-identical pair and nonpair CSVs.

| Structure | Total energy |
| --- | ---: |
| 1RNA | -28.606109849613684 |
| 1EHZ | -133.66638013285746 |

The test also passed for CRLF manifests, a blank row, spaces in output paths,
and an existing empty output directory. Empty/malformed manifests, relative
paths, missing input files/parents, duplicate destinations and nonempty outputs
were rejected before energy-table loading or job output creation.

`tests/test_output_directory.py --docker briq:batch` passed the existing
single-structure five-file export, reference totals, reconstruction, argument
validation and overwrite-protection checks for both fixtures.

## Additional variant sequence

A local check scored the saved 1RNA A:1 variants A/C/G/U, each repeated four
times in a single 16-job process. Every total, pair and nonpair component matched
the previously saved single-run exports within `1e-9`. Observed container memory
after initialization was approximately 2.74-2.80 GiB over this run. This is a short
serial-run check; it is not a memory guarantee for arbitrarily large batches.

## Small timing comparison

The same three jobs took **129.570 s** as separate runs of the previous image,
and **33.741 s** as one batch: **3.84x** faster in this measurement. Times include
Docker startup, initialization, scoring and export. They exclude mutant generation.
Filesystem caching and machine load affect this small comparison; it does not
predict the speedup for a full frustration scan.

## Scope

Batch mode currently fixes seed 0, weights 1/1/1 and PDB phosphate ownership.
Potential tables are shared; structure preparation and scoring remain per job.
This checks implementation compatibility, not the physical validity of the
potential, mutation geometry or a frustration index.

## Total-only batch validation (2026-11-04)

`tests/test_batch_total.py --docker briq:batch-total-test` passed, with the optional
`--variants` directory containing saved local 1RNA A:1 A/C/G/U exports. Eight jobs
were compared: 1RNA, 1EHZ, repeated 1RNA, another 1RNA filename containing spaces
and a comma, and four saved mutants. Both full-batch totals and total-only totals
matched prior saved reference energies with absolute tolerance `1e-9`.

The test checked zero-based manifest ordering, repeated inputs, CSV path escaping,
CRLF and blank rows, one table initialization, no new per-variant output files,
and removal of per-job preparation files from a separately mounted `/tmp`.
Invalid manifests, missing inputs, existing outputs, input/output aliases, missing
output parents, relative output paths and incorrect argument counts were rejected.

The eight-job full-export batch took **68.221 s**, versus **56.947 s** for total-only
(**1.20x**, about 16.5% less elapsed time). This is a small sequential measurement
under the current machine load, not a prediction for a complete frustration scan.
Mutant generation is excluded. Removing exports does not remove full-structure
preparation and scoring; those costs remain.

`tests/test_output_directory.py --docker briq:batch-total-test` also passed the
existing five-file mode for 1RNA and 1EHZ, including argument and overwrite checks.
