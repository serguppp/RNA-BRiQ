# Molecule-wide energy component summaries

`BRiQ_Energy --mutate-total` always writes `energy_components.csv` alongside
`contact_variants.csv`. The existing Python `run_briq_mutations(cache, workspace)`
invocation works unchanged with image `briq:memory`. In that wrapper the new file
is at `workspace/briq/energy_components.csv`.

No per-variant PDB files or pair/nonpair CSV exports are created. Components are
cached with the same interactions as the incremental total score. Their
summation does not alter the original total's accumulation order. The full
reference pair/nonpair export remains available as before.

## Rows and columns

The first row is `state=reference`, `variant_index=-1`, and has empty contact/base
fields. Each subsequent row is `state=variant`, with zero-based `variant_index`
following the row order in `contact_variants.csv`. It repeats `chain_i`,
`residue_i`, `chain_j`, `residue_j`, `contact_type`, `base_i`, and `base_j` for
joining by identity rather than relying on filenames.

Every numeric value is a sum over the **entire molecule**, not just the named
contact. Column meanings match the existing PDB-phosphate-ownership exports:

| Columns | Sum of |
|---|---|
| `E_bb`, `E_bo`, `E_oo`, `E_clash`, `E_pair` | Corresponding column in a full `pairs.csv` export |
| `E_base_ribose`, `E_base_phosphate`, `E_ribose_phosphate`, `E_rotamer`, `E_phosphate_local`, `E_nonpair` | Corresponding column in a full `nonpair.csv` export |
| `E_total` | Original full molecular score |

The following identities hold within floating-point precision:

```text
E_pair = E_bb + E_bo + E_oo + E_clash
E_nonpair = E_base_ribose + E_base_phosphate + E_ribose_phosphate
            + E_rotamer + E_phosphate_local
E_total = E_pair + E_nonpair
```

`E_pair` and `E_nonpair` are subtotals; do not add them again to their constituent
columns. `E_clash` intentionally means the sum of the existing **pair-column**
steric terms. Intra-residue terms remain in the corresponding nonpair columns,
which may themselves include steric contributions. This export preserves the
established bookkeeping rather than redefining clash energy.

To form component differences, subtract the reference row's value from each
variant row's corresponding value. This gives, for example, `delta_E_clash`
without mixing scores from different reference preparations.

Both CSVs are finalized from temporary files at successful completion. The
`mutation_run.json` completion record identifies the sidecar using
`energy_components_file`. Existing output directories are still refused.

## Validation

```bash
python tests/test_mutation_total.py --image briq:memory --verify-full
```

For all 3,728 bundled 1RNA/1EHZ variants this checks unchanged total energies,
matching row identities/order, finite component values, and both subtotal
identities. Native component sums are compared to full CSV exports.
`--verify-full` additionally rebuilds all phosphates and compares all nine
disjoint component sums against the independent full pair/nonpair scorer for
every variant, bypassing the incremental component cache.

Python does not need a signature or return-type change. To show the path in
the existing pipeline, optionally add:

```python
print(f"Energy components: {result.variants_file.with_name('energy_components.csv')}")
```

For a `BriqMutationResult` property, an optional convenience addition is:

```python
@property
def components_file(self) -> Path:
    return self.variants_file.with_name("energy_components.csv")
```
