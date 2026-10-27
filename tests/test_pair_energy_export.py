#!/usr/bin/env python3
"""Integration checks for BRiQ pair-energy export."""

import csv
import math
import pathlib
import re
import shutil
import subprocess
import sys
import tempfile

PAIR_FIELDS = [
    "chain_i",
    "residue_i",
    "chain_j",
    "residue_j",
    "sep",
    "contact_type",
    "E_bb",
    "E_bo",
    "E_oo",
    "E_clash",
    "E_pair",
]
NONPAIR_FIELDS = [
    "chain",
    "residue",
    "E_base_ribose",
    "E_base_phosphate",
    "E_ribose_phosphate",
    "E_rotamer",
    "E_phosphate_local",
    "E_nonpair",
]
TOLERANCE = 1e-9
RUNNER_WORKDIR = None


def invoke(binary, *args):
    command_args = [str(arg) for arg in args]
    if str(binary).startswith("docker:"):
        image = str(binary).split(":", 1)[1]
        container_args = []
        for arg in args:
            if isinstance(arg, pathlib.Path):
                container_args.append(
                    "/work/" + arg.resolve().relative_to(RUNNER_WORKDIR).as_posix()
                )
            else:
                container_args.append(str(arg))
        command = [
            "docker",
            "run",
            "--rm",
            "--mount",
            f"type=bind,source={RUNNER_WORKDIR},target=/work",
            "--workdir",
            "/work",
            image,
            "BRiQ_Energy",
            *container_args,
        ]
    else:
        command = [str(binary), *command_args]
    result = subprocess.run(
        command,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    return result.stdout


def run(binary, *args):
    output = invoke(binary, *args)
    matches = re.findall(r"[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?\s*$", output)
    if not matches:
        raise AssertionError("BRiQ did not report a numeric global energy")
    return float(matches[-1])


def read_csv(path, expected_fields):
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != expected_fields:
            raise AssertionError(f"Unexpected columns in {path}: {reader.fieldnames}")
        return list(reader)


def assert_reconstruction(global_energy, pair_rows, nonpair_rows):
    for row in pair_rows:
        if not all(math.isfinite(float(row[key])) for key in PAIR_FIELDS[6:]):
            raise AssertionError("Non-finite pair component")
    for row in nonpair_rows:
        if not all(math.isfinite(float(row[key])) for key in NONPAIR_FIELDS[2:]):
            raise AssertionError("Non-finite nonpair component")
    reconstructed = sum(float(row["E_pair"]) for row in pair_rows)
    reconstructed += sum(float(row["E_nonpair"]) for row in nonpair_rows)
    if not math.isclose(global_energy, reconstructed, rel_tol=0.0, abs_tol=TOLERANCE):
        raise AssertionError(
            f"Export sums to {reconstructed}, global energy is {global_energy}"
        )


def make_two_chain_fixture(source, destination):
    lines = source.read_text().splitlines()
    atom_lines = [line for line in lines if line.startswith(("ATOM  ", "HETATM"))]
    if not atom_lines:
        raise AssertionError(f"No PDB atoms found in {source}")
    chain_a = []
    chain_b = []
    for line in atom_lines:
        line_a = line[:21] + "A" + line[22:]
        x = float(line[30:38]) + 100.0
        line_b = line[:21] + "B" + line[22:30] + f"{x:8.3f}" + line[38:]
        chain_a.append(line_a)
        chain_b.append(line_b)
    destination.write_text(
        "\n".join(chain_a + ["TER"] + chain_b + ["TER", "END"]) + "\n"
    )


def check_input(binary, pdb_path, workdir, export_option="--pair-energies"):
    workdir.mkdir(parents=True, exist_ok=True)
    input_file = workdir / "input.briq"
    invoke(binary, "--prepare", pdb_path, input_file)
    sequence = next(
        line.split(maxsplit=1)[1]
        for line in input_file.read_text().splitlines()
        if line.startswith("seq ")
    )

    legacy_energy = run(binary, input_file)
    pairs_path = workdir / "pairs.csv"
    nonpair_path = workdir / "nonpair.csv"
    default_energy = run(
        binary,
        input_file,
        export_option,
        pairs_path,
        "--nonpair",
        nonpair_path,
    )
    if not math.isclose(legacy_energy, default_energy, rel_tol=0.0, abs_tol=TOLERANCE):
        raise AssertionError("Default export changed the legacy global score")
    pair_rows = read_csv(pairs_path, PAIR_FIELDS)
    nonpair_rows = read_csv(nonpair_path, NONPAIR_FIELDS)
    expected_pair_count = len(sequence) * (len(sequence) - 1) // 2
    if len(pair_rows) != expected_pair_count:
        raise AssertionError(
            f"Expected {expected_pair_count} unordered pairs, found {len(pair_rows)}"
        )
    pair_keys = [
        (row["chain_i"], row["residue_i"], row["chain_j"], row["residue_j"])
        for row in pair_rows
    ]
    if len(set(pair_keys)) != len(pair_keys):
        raise AssertionError("The export contains duplicate residue pairs")
    if len(nonpair_rows) != len(sequence):
        raise AssertionError("Expected one nonpair row per residue")
    if not any(row["sep"] == "1" for row in pair_rows):
        raise AssertionError("No adjacent-residue pair with sep=1 was exported")
    assert_reconstruction(default_energy, pair_rows, nonpair_rows)

    annotated = ";".join(row["contact_type"] for row in pair_rows)
    if "base_pair" not in annotated or "non_watson_crick_pair" not in annotated:
        raise AssertionError(
            "Both canonical and non-canonical assignments are required"
        )
    if "stacking" not in annotated:
        raise AssertionError(
            "No adjacent stacking contact met the documented criterion"
        )

    weighted_pairs = workdir / "weighted-pairs.csv"
    weighted_nonpair = workdir / "weighted-nonpair.csv"
    weighted_energy = run(
        binary,
        input_file,
        "--pair-energies",
        weighted_pairs,
        "--nonpair",
        weighted_nonpair,
        "--weights",
        "0.5",
        "1.25",
        "0.75",
    )
    assert_reconstruction(
        weighted_energy,
        read_csv(weighted_pairs, PAIR_FIELDS),
        read_csv(weighted_nonpair, NONPAIR_FIELDS),
    )
    if math.isclose(default_energy, weighted_energy, rel_tol=0.0, abs_tol=TOLERANCE):
        raise AssertionError("The second weight configuration did not affect the score")
    return pair_rows


def main():
    global RUNNER_WORKDIR
    if len(sys.argv) == 4 and sys.argv[1] == "--docker":
        binary = "docker:" + sys.argv[2]
        pdb_path = pathlib.Path(sys.argv[3]).resolve()
    elif len(sys.argv) == 3:
        binary = pathlib.Path(sys.argv[1]).resolve()
        pdb_path = pathlib.Path(sys.argv[2]).resolve()
    else:
        raise SystemExit(
            "Usage: test_pair_energy_export.py BRiQ_Energy path/to/1EHZ.pdb\n"
            "   or: test_pair_energy_export.py --docker IMAGE path/to/1EHZ.pdb"
        )
    if (
        not str(binary).startswith("docker:") and not binary.is_file()
    ) or not pdb_path.is_file():
        raise SystemExit("BRiQ_Energy executable and PDB input must exist")

    with tempfile.TemporaryDirectory(prefix="briq-energy-export-") as temp:
        workdir = pathlib.Path(temp)
        RUNNER_WORKDIR = workdir.resolve()
        staged_pdb = workdir / "source.pdb"
        shutil.copyfile(pdb_path, staged_pdb)
        check_input(binary, staged_pdb, workdir)
        two_chain_pdb = workdir / "two-chain.pdb"
        make_two_chain_fixture(staged_pdb, two_chain_pdb)
        rows = check_input(binary, two_chain_pdb, workdir / "two-chain")
        chains = {row["chain_i"] for row in rows} | {row["chain_j"] for row in rows}
        if not {"A", "B"}.issubset(chains):
            raise AssertionError(f"Multi-chain IDs were not retained: {chains}")

    print("BRiQ pair-energy export checks passed")


if __name__ == "__main__":
    main()
