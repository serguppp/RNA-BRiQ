# Scoring fixtures

These PDB files were copied unchanged from the existing RNA Frustratometer
project's `tests/data/`, preserving the inputs used for the original results.
They are not newly downloaded structures or relaxed models.

| Fixture | Residues | Chains | Explicit input break | Original score |
| --- | --- | --- | --- | --- |
| 1RNA | 28 | A/B, 14 each | `break 13` (zero based) | -28.606109849613684 |
| 1EHZ | 76 | A | none | -133.66638013285746 |

The `.briq` files retain the existing sequence and structural annotations.
Only their `pdb` line was made relative. The regression driver rewrites that
line to a container-visible absolute path and writes LF line endings for the
unchanged original BRiQ parser. It records PDB SHA-256 hashes in its report.

1EHZ includes modified nucleotides that original BRiQ standardizes. Its PDB
derived sugar rotamers have zero intrinsic rotamer energy by construction.
See `../../docs/pair-energy-validation.md` for the interpretation of the score,
outgoing phosphates, local terms and contact annotations.
