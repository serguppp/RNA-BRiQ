#!/usr/bin/env python3
"""Validate the supplied RNAComposer PDB, native CIF and Gemmi-converted CIF.

Run with the RNA Frustratometer Python environment (Gemmi is required):
  python tests/test_rnacomposer_inputs.py --data-dir PATH --output NEW_DIRECTORY
The sibling independent-reference helpers and existing Docker images are reused.
This script changes no input structures or scoring implementation.
"""

import argparse
import csv
import hashlib
import json
import math
import runpy
import shutil
import subprocess
import tempfile
from pathlib import Path
from types import SimpleNamespace

import gemmi


def atom_inventory(structure):
    atoms = {}
    for model_index, model in enumerate(structure):
        for chain in model:
            for residue in chain:
                for atom in residue:
                    key = (
                        model_index,
                        chain.name,
                        str(residue.seqid),
                        residue.name,
                        atom.name,
                        atom.altloc,
                    )
                    assert key not in atoms, f"Duplicate atom identity: {key}"
                    atoms[key] = (
                        atom.element.name,
                        (atom.pos.x, atom.pos.y, atom.pos.z),
                    )
    return atoms


def compare_structures(first, second):
    a, b = atom_inventory(first), atom_inventory(second)
    assert a.keys() == b.keys(), "Conversion changed atom/residue/model identities"
    assert all(a[key][0] == b[key][0] for key in a), "Conversion changed elements"
    difference = max(math.dist(a[key][1], b[key][1]) for key in a)
    assert difference == 0, f"Coordinates differ by {difference} A"
    return difference


def summarize_csv(path, numeric):
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    return {
        name: {
            "sum": math.fsum(float(row[name]) for row in rows),
            "minimum": min(float(row[name]) for row in rows),
            "maximum": max(float(row[name]) for row in rows),
            "nonzero_count": sum(float(row[name]) != 0 for row in rows),
        }
        for name in numeric
    }


def native_cif_probe(root, source, image):
    destination = root / "cif_native"
    destination.mkdir()
    shutil.copyfile(source, destination / "input.cif")
    command = [
        "docker",
        "run",
        "--rm",
        "--mount",
        f"type=bind,source={root},target=/work",
        image,
        "BRiQ_Energy",
        "--prepare",
        "/work/cif_native/input.cif",
        "/work/cif_native/input.briq",
    ]
    result = subprocess.run(
        command, capture_output=True, text=True, check=False, timeout=60
    )
    (destination / "stdout.log").write_text(result.stdout)
    (destination / "stderr.log").write_text(result.stderr)
    (destination / "command.json").write_text(json.dumps(command, indent=2))
    prepared = destination / "input.briq"
    has_sequence = prepared.is_file() and any(
        line.startswith("seq ") and len(line.split()) == 2
        for line in prepared.read_text().splitlines()
    )
    # This is an observed input-format limitation, not a successful energy score.
    assert result.returncode != 0 or not has_sequence, (
        "Native CIF preparation now works; extend this regression to score it"
    )
    return {
        "returncode": result.returncode,
        "valid_prepared_sequence": has_sequence,
        "status": "conversion_required_for_this_CIF",
        "stdout": result.stdout,
        "stderr": result.stderr,
        "ATOM_records_matching_original_fixed_prefix": sum(
            line[:6] == "ATOM  " for line in source.read_text().splitlines()
        ),
    }


def check_scored_case(label, case, root, runner, validation, image, reference):
    pdb = case / "input.pdb"
    config = case / "input.briq"
    runner.run(image, "BRiQ_Energy", "--prepare", pdb, config)
    sequence = next(
        line.split()[1]
        for line in config.read_text().splitlines()
        if line.startswith("seq ")
    )
    ids = validation.pdb_ids(pdb)
    assert len(sequence) == len(ids) == 46
    raw = case / "original.tsv"
    runner.run(reference, "BRiQ_Oracle", config, raw)
    oracle = validation.read_rows(raw, "\t")
    original_energy = float(
        next(row["energy"] for row in oracle if row["kind"] == "global")
    )
    pair_path, local_path, components = (
        case / name for name in ["pairs.csv", "nonpair.csv", "components.tsv"]
    )
    result = runner.run(
        image,
        "BRiQ_Energy",
        config,
        "--pair-energies",
        pair_path,
        "--nonpair",
        local_path,
        "--components",
        components,
        "--seed",
        "0",
    )
    (case / "run.stdout.log").write_text(result.stdout)
    (case / "energy.log").write_text(result.stderr)
    validation.close(
        validation.last_energy(result), original_energy, label + " global vs original"
    )
    pairs, report = validation.verify_export(pair_path, local_path, oracle, ids, "pdb")
    report["directional_components"] = validation.compare_raw_components(
        components, oracle
    )
    report["original_global"] = original_energy
    report["exported_global"] = validation.last_energy(result)
    report["sequence"] = sequence
    report["input_config"] = config.read_text().splitlines()
    report["pair_columns"] = summarize_csv(pair_path, validation.PAIR_NUMERIC)
    report["nonpair_columns"] = summarize_csv(local_path, validation.LOCAL_NUMERIC)
    report["representatives"] = validation.representatives(pairs, pdb)
    report["all_ids_unique_complete_and_finite"] = True
    report["input_pdb_sha256"] = hashlib.sha256(pdb.read_bytes()).hexdigest()
    assert report["pair_count"] == 1035 and report["nonpair_count"] == 46
    print(
        f"PASS {label}: E={original_energy:.17g}, pairs=1035, max local error={report['max_pair_component_error']:.3g}",
        flush=True,
    )
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--docker", default="briq:validated")
    parser.add_argument("--reference-image", default="briq:validation-original")
    args = parser.parse_args()
    data = args.data_dir.resolve()
    root = (
        args.output.resolve()
        if args.output
        else Path(tempfile.mkdtemp(prefix="briq-rnacomposer-"))
    )
    root.mkdir(parents=True, exist_ok=True)
    assert not (root / "pdb").exists(), (
        "Choose a new directory to preserve previous outputs"
    )
    validation = SimpleNamespace(
        **runpy.run_path(
            str(Path(__file__).with_name("test_pair_energy_validation.py"))
        )
    )
    sources = {ext: data / f"rnacomposer_example1.{ext}" for ext in ["pdb", "cif"]}
    structures = {ext: gemmi.read_structure(str(path)) for ext, path in sources.items()}
    difference = compare_structures(structures["pdb"], structures["cif"])
    inventory = atom_inventory(structures["pdb"])
    report = {
        "gemmi_version": gemmi.__version__,
        "reference_commit": validation.REFERENCE_COMMIT,
        "image_ids": {
            image: subprocess.check_output(
                ["docker", "image", "inspect", image, "--format", "{{.Id}}"], text=True
            ).strip()
            for image in [args.docker, args.reference_image]
        },
        "source_sha256": {
            ext: hashlib.sha256(path.read_bytes()).hexdigest()
            for ext, path in sources.items()
        },
        "input_comparison": {
            "models": len(structures["pdb"]),
            "chains": [chain.name for chain in structures["pdb"][0]],
            "residues": sum(len(chain) for chain in structures["pdb"][0]),
            "atoms": len(inventory),
            "heavy_atoms": sum(
                value[0] not in {"H", "D"} for value in inventory.values()
            ),
            "identical_atom_ids_and_elements": True,
            "maximum_coordinate_difference": difference,
        },
        "native_cif": native_cif_probe(root, sources["cif"], args.docker),
        "cases": {},
    }
    (root / "report.json").write_text(json.dumps(report, indent=2))
    print("Artifacts:", root, flush=True)
    print("Native CIF:", report["native_cif"]["status"], flush=True)
    for label in ["pdb", "cif_to_pdb"]:
        case = root / label
        case.mkdir()
        if label == "pdb":
            shutil.copyfile(sources["pdb"], case / "input.pdb")
        else:
            structures["cif"].write_pdb(str(case / "input.pdb"))
            compare_structures(
                structures["cif"], gemmi.read_structure(str(case / "input.pdb"))
            )
        report["cases"][label] = check_scored_case(
            label,
            case,
            root,
            validation.Runner(root),
            validation,
            args.docker,
            args.reference_image,
        )
        # Keep command logs separate for each case; Runner starts from one.
        for log in root.glob("command-*.log"):
            shutil.move(str(log), str(case / log.name))
        (root / "report.json").write_text(json.dumps(report, indent=2))
    for filename in ["pairs.csv", "nonpair.csv", "components.tsv"]:
        assert (root / "pdb" / filename).read_bytes() == (
            root / "cif_to_pdb" / filename
        ).read_bytes(), filename
    assert [
        line
        for line in (root / "pdb/input.briq").read_text().splitlines()
        if not line.startswith("pdb ")
    ] == [
        line
        for line in (root / "cif_to_pdb/input.briq").read_text().splitlines()
        if not line.startswith("pdb ")
    ]
    runner = validation.Runner(root)
    repeat = root / "repeat.csv"
    result = runner.run(
        args.docker,
        "BRiQ_Energy",
        root / "pdb/input.briq",
        "--pair-energies",
        repeat,
        "--seed",
        "0",
    )
    validation.close(
        validation.last_energy(result),
        report["cases"]["pdb"]["original_global"],
        "Repeat score",
    )
    assert repeat.read_bytes() == (root / "pdb/pairs.csv").read_bytes()
    assert (
        Path(str(repeat) + ".nonpair.csv").read_bytes()
        == (root / "pdb/nonpair.csv").read_bytes()
    )
    report["converted_CIF_and_original_PDB_exports_byte_identical"] = True
    report["repeat_byte_identical"] = True
    report["status"] = "passed_with_CIF_conversion_required"
    (root / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(
        "PASS: all local/original comparisons, conversion equivalence and repeatability; native CIF requires conversion.",
        flush=True,
    )


if __name__ == "__main__":
    main()
