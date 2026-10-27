# BRiQ pair-energy validation

Scope: validate the pair-energy exporter against original BRiQ. No new potential, mutation generator or frustration implementation. Numerical results and representative contacts are in `validation-results.json` alongside this report.

**Outcome (2026-10-23): passed for both supplied fixtures, including every local pair component, ordered helper contribution and seeded weighted repeat.** The original BRiQ score is preserved. This is validation of faithful decomposition under the documented ownership convention, not independent validation of BRiQ's scientific accuracy.

## Reference

The implementation is based on commit `2684f6e6b2526bfbe4ca9ac31773f496384d0b4d`
in [Tomasz Zok's fork](https://github.com/tzok/RNA-BRiQ), including the energy
executable from [upstream PR 4](https://github.com/Jian-Zhan/RNA-BRiQ/pull/4).
The archived reference sources were compared with the upstream implementation
before validation. Results below document the validation performed on 2026-10-23.

## Corrections

1. **Phosphate ownership:** node `j` stores the outgoing phosphate between riboses `j` and `j+1`. `BRTreeInfo::printPDB()` writes its P/O5'/OP1/OP2 atoms to PDB residue `j+1`. The previous exporter labeled that phosphate as residue `j`. The total was correct but physical pair identities were wrong. Base-phosphate and ribose-phosphate now map to `(i,j+1)`; phosphate-phosphate maps to `(i+1,j+1)`, only for existing links. Coincident physical endpoints become same-residue nonpair terms. Original helper calls and total accumulation are unchanged.
2. **Compatibility:** CSV columns are unchanged. Default `--pair-ownership pdb` uses PDB atom ownership. `--pair-ownership node` reproduces the previous grouping. The mode is reported on stderr before the final energy line. Existing consumers requiring old pair values must request `node` explicitly.
3. **Precision:** component TSV now uses 17 significant digits; six-digit stream output cannot support `1e-9` comparisons.
4. **Chains:** `--prepare` preserves geometric break detection and also inserts breaks at changes of PDB chain ID.
5. **Deterministic weighting:** `--seed N` fixes original loop-closure RNG choices. The legacy default remains time-based. These choices matter when BREAK and CONNECT differ; they are distinct from physical input `break` positions.
6. **Checks:** fixed the old test invocation, tightened finite/numerical checks, added independent reference compilation and all-pair regression tests. Non-finite global energy is rejected by the executable.

No files under `para/`, `forcefield/`, `geometry/`, or `model/` were modified. Parameter data, weights, mathematical formulas and original coordinate preparation were preserved.

## Complete scoring path

Let A/B be original BRiQ nodes, P(B) the outgoing phosphate of B, q the outer CLASH weight and wO=1.2. Distances are in angstroms. Coordinates are those prepared by the original constructor, including standardized bases and reconstructed phosphates.

| Raw TSV term | Meaning and atoms | Weight and zero conditions | Pair CSV |
| --- | --- | --- | --- |
| `base_base` | Statistical 6D base-frame interaction, including pairing and stacking | Internal separation weights `.4/.8/1`; adjacent pyrimidine-pyrimidine factor `.5`. Helper requires anchor distance squared `<200` and minimum base-atom distance `<4.5`. Table also rejects minimum distance `<1.5` or `>=5`, frame distance `>=15`, missing entries. No q. | `E_bb(A,B)`, once for A<B |
| `base_clash` | Base-base steric lookup, sampled in both frames within one helper call | q; self, distant anchors (`>=16`), outside grids or beyond clash thresholds give zero | `E_clash(A,B)`, once; do not add the reversed helper again |
| `ribose_clash` | Combined 8-by-8 sugar repulsion **and** O2'-O2' H-bond | q multiplies both parts. Outer sugar atom-2 distance `<9`; H-bond requires `abs(sep)>=2` and O2'-O2' `<4`. Atomic tables reject distances `>=4`; zero table entries can eliminate either part. | H-bond to `E_oo(A,B)`, repulsion to `E_clash(A,B)`, once |
| `phosphate_clash` | Outgoing O5'/OP1/OP2 against O5'/OP1/OP2 | q; both links required; P-P `<=6`; helper loops over atom distances `<5`, with shorter cutoffs in atomic tables | `E_clash(P(A),P(B))`, both endpoints shifted |
| `base_ribose` | Directional base against O2'/O3'/O4' tables, plus base against 8 sugar atoms | `q*(wO*oxygen+repulsion)`; `sep!=0`, anchor-to-sugar atom-2 `<12`. Oxygen has grid, separation/type and angular gates. Atomic repulsion cutoff 4. | Both original ordered calls aggregated: oxygen to `E_bo`, repulsion to `E_clash` |
| `base_phosphate` | Directional base against O5'/OP1/OP2 tables and repulsion | `q*(wO*oxygen+repulsion)`; existing outgoing link, `sep!=0`; same oxygen lookup gates; atomic repulsion cutoff 4 | Physical `(A,P(B))`: oxygen to `E_bo`, repulsion to `E_clash`; coincident endpoints to nonpair |
| `ribose_phosphate` | 8 sugar atoms against O5'/OP1/OP2, plus O2'-nearest-OP H-bond | `q*(Hbond+repulsion)`; absent phosphate, `sep==0` or `-1` returns zero. H-bond only `sep>2`, nearest OP `<4`; atom repulsion cutoff 4. | Physical `(A,P(B))`: H-bond to `E_oo`, repulsion to `E_clash`; coincident endpoints to nonpair |
| `rotamer` | Intrinsic sugar/rotamer preference | Library energies already scaled by `.7`; complete PDB-backbone rotamers explicitly receive zero. No q. | Nonpair `E_rotamer`, anchored at A |
| `phosphate_local` | Backbone link geometry from **two adjacent riboses**: bond/angle penalties and statistical torsions | `pho.e` includes `wtPho=.7`, `wtDihed=1.2`, `kBond=5`, `kAng=.1`. Outer BREAK/CONNECT chosen by `findMidPoint()` flags. Missing outgoing link gives default zero. No q. | Nonpair `E_phosphate_local`, anchored upstream; link-local, **not physically one-body** |

These are all terms in `totalEnergy()`. A separate `BaseStackingEnergyTable` and `getPhoPhoHbondEnergy()` exist but are not called on this path; adding them would change the score. Atomic clash helpers remain original table lookups. BRiQ score units are preserved; no kcal/mol conversion is asserted.

The atom tables are generated with original `vdwEnergy()`, internal `wtClash=1` and `lamdaClash=2.8`: zero beyond the atom-pair clash distance, quartic inside it, and a linear continuation for deep overlaps. The base-frame clash helper includes its original half factors. The H-bond tables also retain original loader shifts/rescaling (`O2O2: +2.5`, `O2OP: +3.5` before `energyRescale`). These internal operations are already included in helper returns and must not be weighted a second time.

The raw component TSV retains original one-based **node** indices and ordered helper values. `ribose_clash` therefore remains a combined helper term there. In PDB ownership mode CSV `sep` describes the reported physical endpoint residues; a phosphate contribution may have been computed with a different node separation. Do not recompute phosphate components from CSV `sep` alone.

`E_pair = E_bb + E_bo + E_oo + E_clash`. Same-physical-residue terms, rotamers and phosphate-link local energies are separately reported in the nonpair file. One nonpair row per residue is bookkeeping, not a claim that all terms depend on one residue only.

## Sequence, breaks and contacts

Input `break k` removes the outgoing covalent link after zero-based node k. For 1RNA `break 13` separates two 14-residue chains; neither chain end has an outgoing phosphate. Separation is 0 for self, +/-1 for connected neighbors, +/-2 for two connected steps and **3 in either direction** otherwise, including across breaks. Index adjacency alone does not establish covalent adjacency.

`seq` supplies the input length and must agree with the PDB; actual base types come from parsed, standardized PDB bases. Changing only `seq` letters is not mutation. `sec` and `nwc` build the sampling tree and select loop flags. They do not directly add WC/non-WC energy bonuses, but can affect which local phosphate terms receive unequal BREAK/CONNECT weights. Annotations cannot generally be described as irrelevant to scoring.

`contact_type` is independent geometric/AssignRNASS annotation, not an energy threshold. Its stacking rule currently detects only connected adjacent residues. `unannotated` does not certify absence of contact. Contacts may have zero tabulated energy, while backbone terms may contribute without a base-base contact. A frustration contact mask must be defined independently of `E_pair != 0`.

Concrete examples are 1RNA A:4-B:26 and 1EHZ A:1-A:71: both satisfy the exported geometric `base_contact` rule yet have `E_pair=0`. There are respectively 1 and 33 geometrically annotated zero-energy pairs in these fixtures.

## Measured zero components in 1RNA

The original oracle gives `E_bo=E_oo=0`, before CSV aggregation; no cancellation or artificial replacement is involved. The actual base-oxygen weight is **1.2** and outer CLASH is **1**, so the terms are active.

Of 2640 original base-oxygen atom queries reaching the relevant helpers, 180 are excluded by separation/type rules, 2236 lie outside the table grid, and 224 select a zero table entry. None produces a nonzero value; no angular rejection is needed to explain this fixture.

For O2'-O2', the 26 connected neighbor pairs are excluded from H-bond evaluation. All other 352 pairs fail the outer C3'-C3' `<9 A` gate; direct measurement of the complete PDB sugars confirms this. For O2'-OP, 122 ordered combinations with existing phosphates are excluded by separation, and all remaining 606 have nearest OP distance `>=4 A`. The code computes these terms in `totalEnergy()`; this structure's geometry and original cutoffs give zero.

1EHZ is a positive control: 34 original base-ribose and 16 base-phosphate ordered helper calls have nonzero oxygen contributions. One ribose-ribose and two ribose-phosphate helper calls have nonzero H-bond contributions. After physical ownership aggregation there are 36 pairs with nonzero `E_bo` and 3 with nonzero `E_oo`.

## Independent regression method

`tests/energy_oracle.cpp` is compiled against a fresh archive of the fixed unmodified reference commit, including its original header-only helper functions. It never links to the modified scorer. It calls every original unordered base/sugar/phosphate helper and every ordered base-sugar/base-phosphate/sugar-phosphate helper. Oxygen and H-bond contributions are queried directly from underlying original tables, independently of Copilot's optional output pointers. Python independently aggregates original directions and maps atom ownership, then checks every exported component with absolute tolerance `1e-9`.

The driver checks fixed original global baselines, reconstruction, unique complete pair sets, metadata, insertion codes and repeated residue numbers across chains, finite values, raw TSV directions, both ownership modes, weighted scores with seed 0, byte-identical repeats and the original one-argument invocation. Fixtures were copied from existing project data. Generated Linux inputs use LF: original `InputParser` does not strip CR from Windows CRLF lines.

After building the candidate image with the BRiQ parameter data, run from the repository root:

```powershell
python tests/test_pair_energy_validation.py --docker briq:validated --build-reference
python tests/test_pair_energy_export.py --docker briq:validated tests/fixtures/1EHZ.pdb
```

`--build-reference` compiles image `briq:validation-original` from the fixed original commit; the candidate supplies the parameter files and compiler tools, not the reference code. The existing parent Dockerfile provides these dependencies. The driver prints its artifact directory; `--output PATH` selects a new persistent directory. All logs, original helper values, diagnostics and CSVs remain there. The existing `briq:tzok` image was not overwritten.

## Results

All comparisons used absolute tolerance `1e-9`. Original and exported global values were identical in both cases.

| Structure | Original/exported global | Pair + nonpair sum | Pairs / nonpair rows | Maximum pair-component difference |
| --- | --- | --- | --- | --- |
| 1RNA | -28.606109849613684 | -28.606109849613702 | 378 / 28 | 0 |
| 1EHZ | -133.66638013285746 | -133.66638013285734 | 2850 / 76 | 4.440892098500626e-16 |

Nonpair component differences were zero. All 3920 raw component records for 1RNA and 28880 for 1EHZ matched the original ordered/unordered helper returns exactly. Reconstruction errors were respectively `1.7763568394002505e-14` and `1.1368683772161603e-13`. No NaN/Inf, duplicate/reversed pair, missing ID or unexpected record count was found.

With outer weights BREAK=.5, CONNECT=1.25, CLASH=.75 and seed 0, reconstructed energies were `-9.51211987333464` (1RNA) and `-118.78543566813352` (1EHZ), agreeing with separately evaluated original weighted scores within tolerance. Default and seeded weighted reruns produced byte-identical pair and nonpair CSVs. Original one-argument scoring also matched. The legacy integration script passed its original and synthetic two-chain 1EHZ cases. Ruff lint/format and Git whitespace checks passed.

Representative local energies below are **physical PDB ownership**, each checked against original functions independently of CSV sums. Complete component values are in the JSON report.

| Structure / class | Residue pair | E_bb | E_bo | E_oo | E_clash | E_pair |
| --- | --- | --- | --- | --- | --- | --- |
| 1RNA adjacent stacking | A:1-A:2 | -0.87947599097674267 | 0 | 0 | 0 | -0.87947599097674267 |
| 1RNA WC | A:1-B:28 | -2.9151526362090809 | 0 | 0 | 0 | -2.9151526362090809 |
| 1RNA positive clash | A:5-A:6 | -2.2162197916796518 | 0 | 0 | 0.00029253475788455194 | -2.2159272569217672 |
| 1EHZ adjacent stacking | A:1-A:2 | -2.7862799739010358 | 0 | 0 | 0 | -2.7862799739010358 |
| 1EHZ WC | A:1-A:72 | -6.8111034314913699 | 0 | 0 | 0 | -6.8111034314913699 |
| 1EHZ non-WC | A:8-A:14 | -3.2694080333241438 | -0.0204 | 0 | 0 | -3.2898080333241437 |
| 1EHZ backbone/clash | A:7-A:49 | 0 | -0.2952 | 0.4700000000000002 | 0.000026251826009960548 | 0.17482625182601019 |

The non-contact controls A:1-A:6 have zero in every component for both structures; independent minimum input heavy-atom distances are 14.8333 A (1RNA) and 13.2242 A (1EHZ). 1RNA has no non-WC assignment in the supplied input; the noncanonical representative is taken from 1EHZ rather than invented.

Phosphate remapping changes **44** pair-component records for 1EHZ and none for 1RNA (whose phosphate pair contributions are zero). For example, the previous node grouping gave A:7-A:48 `E_pair=0.47002625182601016`; the correct physical grouping gives that pair zero and attributes the outgoing phosphate contribution to A:7-A:49. These changes preserve the global sum.

The new frozen local checks were also executed against the **preserved Copilot image** `briq:validation-copilot`. Its global energy still matched, but the local regression failed as intended: A:8-A:14 was `-3.2694080333241438` instead of `-3.2898080333241437`, and A:7-A:48 had an erroneously attributed `E_oo=0.4700000000000002`. The corrected exporter passes both checks.

A second regression retains covalent coordinates but changes chain IDs in the middle of 1RNA. `--prepare` now emits the required extra break after node 5, as well as break 13. The preserved earlier-exporter image failed this check; the corrected image passed. Renumbering residues independently per chain and adding insertion code `2A` also preserved energies and correct IDs.

The complete final run and the additional artifact audit passed. A preliminary run using CRLF test input was rejected because the unchanged original parser retained CR; test generation was corrected to LF. A later redundant preliminary run was deliberately stopped after its default local comparisons passed, to free memory for the complete final run. Neither interrupted run is used as evidence of a completed suite.

## Limitations and next stage

Additional RNAComposer testing (2026-10-23) validated the supplied 46-residue PDB and a Gemmi conversion of its atom-identical CIF: both score `37.066634843230638`, with maximum pair-component error `1.0408340855860843e-16` and identical exports. The supplied CIF **fails direct native parsing**, including an observed SIGSEGV; use Gemmi conversion for that input. See `tests/test_rnacomposer_inputs.py` to reproduce these checks with externally supplied RNAComposer fixtures. No scoring-engine changes were made for these additional tests.

Agreement is with the specified original implementation and parameter set, not proof of experimental accuracy or exhaustive geometric coverage. Original code standardizes modified nucleotides, may choose library rotamers for incomplete backbones, and reconstructs phosphate coordinates. PDB-based geometric labels need not exactly match the prepared scoring geometry.

Aggregate CSV loses directional subterms and their original node separation; use raw TSV/oracle outputs for audits. A later mutation calculation must preserve original helper dependencies and explicitly decide whether its target is a pair energy or a neighborhood including nonpair link terms.

Default RNG is time-based. Equal BREAK/CONNECT weights remove this particular weight-selection dependence for these fixtures; use `--seed` for reproducible nondefault weighting. This does not validate every possible upstream stochastic preparation path.

Some upstream invalid-input paths return `exit(0)`. Tests check required files and a finite energy instead of relying solely on the return code. General parser hardening and unsupported/incomplete input cases remain outside this validation.

An inherited scientific-code detail remains uncertain: `AtomicEnergyTable::getPhoPhoEnergy()` indexes the common atomic tables with `typeA*36+typeB`, whereas other atom helpers use `typeA*47+typeB` and the loader creates 2209 (=47*47) rows. This is present in the unmodified reference. It was **not corrected or declared intentional**: changing it would change the original potential/score, contrary to this task. The exporter faithfully reproduces its returned phosphate-clash energies. Confirming the intended phosphate atom mapping requires a separate upstream investigation.

Before a later 16-variant A/C/G/U experiment, define which base coordinates, sugar rotamers, reconstructed phosphates and structural annotations remain fixed or are rebuilt. Mutation generation, relaxation, reference ensembles and frustration statistics are not implemented or validated here.

The exporter is ready as the validated scoring/decomposition baseline for designing that next experiment on these supported inputs. A universally validated mutational protocol or an exhaustive RNA parser is not claimed. Keep the physical ownership convention fixed, preserve an independent contact mask, and retain nonpair terms explicitly when comparing variants.
