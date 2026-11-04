# RNAfrustratometer Python integration

The BRiQ implementation is available in image `briq:memory`. These are proposed
changes for the separate Python RNAfrustratometer application; they have not
been applied to that application's source files.

## Add to briq.py

Merge the imports below with the existing imports, then add the dataclass and
function. `MutationCache` is imported from the RNAfrustratometer package.
The function exports one prepared reference and one geometry cache, invokes
BRiQ once, and validates every returned contact's 16 total energies.

```python
"""Code to copy into briq.py; application source files are not modified."""
import csv
import json
import math
import subprocess
from dataclasses import dataclass
from itertools import product
from pathlib import Path

from rna_frustratometer.mutation import MutationCache


@dataclass(frozen=True)
class BriqMutationResult:
    total_energy: float
    pairs_file: Path
    variants_file: Path
    run_file: Path


def run_briq_mutations(
    cache: MutationCache,
    workspace: Path,
    *,
    image: str = "briq:memory",
) -> BriqMutationResult:
    workspace = workspace.resolve()
    workspace.mkdir(parents=True, exist_ok=True)
    reference = workspace / "reference.pdb"
    geometry = workspace / "geometry.tsv"
    output = workspace / "briq"

    if reference.exists() or geometry.exists() or output.exists():
        raise FileExistsError("Use a new analysis output directory.")

    cache.reference.write_pdb(str(reference))
    with geometry.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["chain", "residue", "base", "atom", "x", "y", "z"])
        for (chain, identifier, base), residue in sorted(cache.variants.items()):
            if residue.name != base:
                raise ValueError("Unexpected cached residue name.")
            for atom in residue:
                if atom.element.name in ("H", "D"):
                    continue
                if atom.altloc != "\x00":
                    raise ValueError("Select one atom conformation first.")
                if not all(math.isfinite(v) for v in (atom.pos.x, atom.pos.y, atom.pos.z)):
                    raise ValueError("Non-finite atom coordinate.")
                writer.writerow([
                    chain, identifier, base, atom.name,
                    f"{atom.pos.x:.3f}", f"{atom.pos.y:.3f}", f"{atom.pos.z:.3f}",
                ])

    with (workspace / "memory.log").open("w", encoding="utf-8") as log:
        subprocess.run(
            [
                "docker", "run", "--rm", "--mount",
                f"type=bind,source={workspace},target=/work",
                image, "BRiQ_Energy", "--mutate-total",
                "/work/reference.pdb", "/work/geometry.tsv", "/work/briq",
            ],
            check=True,
            stdout=log,
            stderr=subprocess.STDOUT,
        )

    run_file = output / "mutation_run.json"
    metadata = json.loads(run_file.read_text(encoding="utf-8"))
    native = float(metadata["native_energy"])
    if not math.isfinite(native):
        raise ValueError("Invalid native energy.")
    variants_file = output / "contact_variants.csv"
    with variants_file.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        expected = [
            "chain_i", "residue_i", "chain_j", "residue_j", "contact_type",
            "base_i", "base_j", "native_E_total", "E_total", "delta_E_total",
        ]
        if reader.fieldnames != expected:
            raise ValueError("Invalid mutation CSV header.")
        rows = list(reader)

    groups: dict[tuple[str, str, str, str], list[tuple[str, str]]] = {}
    for row in rows:
        key = (row["chain_i"], row["residue_i"], row["chain_j"], row["residue_j"])
        groups.setdefault(key, []).append((row["base_i"], row["base_j"]))
        energy = float(row["E_total"])
        reference_energy = float(row["native_E_total"])
        delta = float(row["delta_E_total"])
        if not all(math.isfinite(v) for v in (energy, reference_energy, delta)):
            raise ValueError("Invalid mutant energy.")
        if abs(reference_energy - native) > 1e-9 or abs(delta - (energy-native)) > 1e-9:
            raise ValueError("Inconsistent reference or energy difference.")
    combinations = list(product("ACGU", repeat=2))
    if any(variants != combinations for variants in groups.values()):
        raise ValueError("Expected 16 ordered variants per contact.")
    if len(rows) != metadata["variants"] or len(groups) != metadata["contacts"]:
        raise ValueError("Mutation results do not match run metadata.")

    return BriqMutationResult(
        total_energy=native,
        pairs_file=output / "pairs.csv",
        variants_file=variants_file,
        run_file=run_file,
    )
```

## Replace RNAfrustratometer.run() in pipeline.py

Merge these imports:

```python
from shutil import copyfile
from .briq import run_briq_mutations
```

Then replace the method:

```python
    def run(self) -> AnalysisResult:
        if not self.input_dir.is_file():
            raise FileNotFoundError(self.input_dir)
        if self.input_dir.suffix.lower() != ".pdb":
            raise ValueError("Input must be a PDB file.")

        self.load_input_structure()
        if self.input_structure is None:
            raise RuntimeError("Input structure not loaded.")
        if len(self.input_structure) != 1:
            raise ValueError("Select exactly one structural model.")

        # Contacts will be discovered by the same C++ process that scores them.
        # Prepare every residue now, so no separate native-scoring call is needed.
        cache = MutationCache(self.input_structure)
        targets = {
            (chain.name, str(residue.seqid))
            for chain in self.input_structure[0]
            for residue in chain
        }
        cache.prepare(targets)
        result = run_briq_mutations(cache, workspace=self.output_dir)

        self.input_briq = BriqResult(
            total_energy=result.total_energy,
            pairs_file=result.pairs_file,
            nonpair_file=result.pairs_file.with_name("nonpair.csv"),
            run_file=result.pairs_file.with_name("run.json"),
        )
        self.load_contacts()
        energies: dict[PairKey, list[float]] = {key: [] for key, _ in self.contacts}

        with result.variants_file.open(newline="", encoding="utf-8") as file:
            for row in csv.DictReader(file):
                key = (
                    row["chain_i"], row["residue_i"],
                    row["chain_j"], row["residue_j"],
                )
                energies[key].append(float(row["E_total"]))

        contact_results = tuple(
            ContactResult(
                chain_i=key[0],
                residue_i=key[1],
                chain_j=key[2],
                residue_j=key[3],
                native_energy=result.total_energy,
                variant_energies=tuple(energies[key]),
            )
            for key, _ in self.contacts
        )
        if any(len(contact.variant_energies) != 16 for contact in contact_results):
            raise ValueError("Expected 16 energies for each contact.")

        summary = self.output_dir / "contact_variants.csv"
        copyfile(result.variants_file, summary)
        print(f"Original energy: {result.total_energy:.6f}")
        print(f"Contacts: {len(contact_results)}")
        print(f"Saved: {summary}")
        return AnalysisResult(
            native_energy=result.total_energy,
            contacts=contact_results,
        )
```

Keep the existing `ContactResult.frustration` formula, `AnalysisResult`,
`load_input_structure`, `read_pairs` and `load_contacts`. The old
`load_input_briq` method is no longer called by this pipeline; no preliminary
`run_briq` or `run_briq_totals` call is needed. Remove the unused
`run_briq_totals` import. `MutationCache.write_pair_mutants` remains available
for debugging, but is not called by this pipeline.

Use a new output directory for every run. Inputs are `reference.pdb` and
`geometry.tsv` in the workspace; BRiQ outputs are in `briq/`, and the pipeline
copies the summary to `contact_variants.csv` at the workspace root.

```python
from pathlib import Path
from rna_frustratometer.pipeline import RNAfrustratometer

result = RNAfrustratometer(
    input_dir=Path("tests/data/1EHZ.pdb"),
    output_dir=Path("results/1EHZ_memory_01"),
).run()
```

The method and helper have been executed in a temporary Python process against
both structures, without editing application source. Full per-variant PDB files
are no longer needed. Geometry is still supplied by the existing ModeRNA cache;
this optimization does not change the scientific frustration definition.
