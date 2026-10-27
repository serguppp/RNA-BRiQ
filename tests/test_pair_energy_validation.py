#!/usr/bin/env python3
"""Compare every exported pair with an oracle built from unmodified BRiQ.

Run from the repository root:
  python tests/test_pair_energy_validation.py --docker briq:validated --build-reference
Only Docker and Python's standard library are required. No frustration is computed.
"""

import argparse
import csv
import hashlib
import io
import itertools
import json
import math
import shutil
import subprocess
import tempfile
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

REFERENCE_COMMIT = "2684f6e6b2526bfbe4ca9ac31773f496384d0b4d"
EXPECTED = {"1RNA": -28.606109849613684, "1EHZ": -133.66638013285746}
# Frozen local references supplement the independently compiled all-pair oracle.
# The last two 1EHZ contacts catch the outgoing-phosphate ownership error.
LOCAL_BASELINES = {
    "1RNA": {
        ("A", "1", "A", "2"): -0.87947599097674267,
        ("A", "1", "B", "28"): -2.9151526362090809,
        ("A", "1", "A", "6"): 0.0,
        ("A", "5", "A", "6"): -2.2159272569217672,
    },
    "1EHZ": {
        ("A", "1", "A", "2"): -2.7862799739010358,
        ("A", "1", "A", "72"): -6.8111034314913699,
        ("A", "8", "A", "14"): -3.2898080333241437,
        ("A", "1", "A", "6"): 0.0,
        ("A", "7", "A", "48"): 0.0,
        ("A", "7", "A", "49"): 0.17482625182601019,
    },
}
TOL = 1e-9
PAIR_NUMERIC = ["E_bb", "E_bo", "E_oo", "E_clash", "E_pair"]
LOCAL_NUMERIC = [
    "E_base_ribose",
    "E_base_phosphate",
    "E_ribose_phosphate",
    "E_rotamer",
    "E_phosphate_local",
    "E_nonpair",
]
RAW_KINDS = {
    "base_base",
    "base_clash",
    "ribose_clash",
    "phosphate_clash",
    "base_ribose",
    "base_phosphate",
    "ribose_phosphate",
    "rotamer",
    "phosphate_local",
}


def close(actual, expected, message):
    if (
        not math.isfinite(actual)
        or not math.isfinite(expected)
        or abs(actual - expected) > TOL
    ):
        raise AssertionError(f"{message}: {actual:.17g} != {expected:.17g}")


def read_rows(path, delimiter=","):
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter=delimiter))


def check_local_baselines(name, pairs):
    indexed = {
        tuple(row[key] for key in ["chain_i", "residue_i", "chain_j", "residue_j"]): row
        for row in pairs
    }
    for key, expected in LOCAL_BASELINES[name].items():
        close(float(indexed[key]["E_pair"]), expected, f"{name} fixed local {key}")


def pdb_ids(path):
    """Independent PDB author IDs, including insertion codes, for RNA residues."""
    groups = {}
    models = 0
    for line in path.read_text().splitlines():
        if line.startswith("MODEL "):
            models += 1
            if models > 1:
                break
        if not line.startswith(("ATOM  ", "HETATM")) or len(line) < 54:
            continue
        key = (line[21:22], line[22:27].strip())
        groups.setdefault(key, set()).add(line[12:16].strip())
    return [key for key, names in groups.items() if "C1'" in names]


def build_reference(repo, work, image, reference):
    context = work / "reference-build"
    context.mkdir()
    archive = subprocess.check_output(
        ["git", "-C", str(repo), "archive", "--format=zip", REFERENCE_COMMIT]
    )
    with zipfile.ZipFile(io.BytesIO(archive)) as zipped:
        zipped.extractall(context / "original")
    shutil.copyfile(repo / "tests/energy_oracle.cpp", context / "energy_oracle.cpp")
    (context / "Dockerfile").write_text(
        f"FROM {image}\nCOPY original /opt/briq-original\nCOPY energy_oracle.cpp /opt/energy_oracle.cpp\n"
        "RUN cmake -S /opt/briq-original -B /opt/briq-original/build >/tmp/reference-config.log 2>&1 "
        "&& cmake --build /opt/briq-original/build --parallel 4 >/tmp/reference-build.log 2>&1 "
        "&& g++ -O2 -std=c++11 -pthread -I/opt/briq-original /opt/energy_oracle.cpp "
        "/opt/briq-original/build/lib/libbriq.a -o /usr/local/bin/BRiQ_Oracle\n"
    )
    with (work / "reference-build.log").open("w") as log:
        subprocess.run(
            ["docker", "build", "-t", reference, str(context)],
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )


class Runner:
    def __init__(self, root):
        self.root = root
        self.number = 0

    def run(self, image, executable, *args):
        converted = [
            "/work/" + arg.relative_to(self.root).as_posix()
            if isinstance(arg, Path)
            else str(arg)
            for arg in args
        ]
        command = [
            "docker",
            "run",
            "--rm",
            "--mount",
            f"type=bind,source={self.root},target=/work",
            "--workdir",
            "/work",
            image,
            executable,
            *converted,
        ]
        self.number += 1
        result = subprocess.run(
            command,
            text=True,
            capture_output=True,
            check=False,
            timeout=300,
        )
        (self.root / f"command-{self.number:03d}.log").write_text(
            "COMMAND "
            + repr(command)
            + "\nSTDOUT\n"
            + result.stdout
            + "\nSTDERR\n"
            + result.stderr
        )
        if result.returncode:
            raise AssertionError(
                f"Command failed ({result.returncode}): {command}\n{result.stderr[-3000:]}\n{result.stdout[-1000:]}"
            )
        return result


def last_energy(result):
    try:
        value = float(result.stderr.strip().splitlines()[-1].split()[-1])
    except (ValueError, IndexError) as exc:
        raise AssertionError("No final numeric energy on stderr") from exc
    if not math.isfinite(value):
        raise AssertionError("Non-finite global energy")
    return value


def check_chain_boundary(pdb, runner, image):
    """A chain-ID boundary must break even geometrically continuous RNA."""
    ids = pdb_ids(pdb)
    reassigned = set(ids[6:14])
    changed = []
    for line in pdb.read_text().splitlines():
        if (
            line.startswith(("ATOM  ", "HETATM"))
            and (line[21:22], line[22:27].strip()) in reassigned
        ):
            line = line[:21] + "C" + line[22:]
        changed.append(line)
    split = pdb.parent / "chain-split.pdb"
    split.write_text("\n".join(changed) + "\n", newline="\n")
    prepared = pdb.parent / "chain-split.briq"
    runner.run(image, "BRiQ_Energy", "--prepare", split, prepared)
    breaks = {
        int(value)
        for line in prepared.read_text().splitlines()
        if line.startswith("break ")
        for value in line.split()[1:]
    }
    assert {5, 13} <= breaks, "A chain-ID boundary was treated as a covalent link"


def aggregate_oracle(rows, size, ownership):
    pairs = {
        key: dict.fromkeys(PAIR_NUMERIC, 0.0)
        for key in itertools.combinations(range(size), 2)
    }
    local = [dict.fromkeys(LOCAL_NUMERIC, 0.0) for _ in range(size)]
    phosphate = {
        int(row["i"]): int(row["j"]) for row in rows if row["kind"] == "phosphate_owner"
    }

    def owner(index):
        return phosphate[index] if ownership == "pdb" else index

    def add(i, j, column, value, local_column=None):
        if i < 0 or j < 0:
            close(value, 0, "Absent phosphate has a nonzero interaction")
        elif i == j:
            if local_column is None:
                close(value, 0, "Unexpected intra-residue interaction")
            else:
                local[i][local_column] += value
        else:
            pairs[tuple(sorted((i, j)))][column] += value

    for row in rows:
        kind = row["kind"]
        if kind not in RAW_KINDS:
            continue
        i, j = int(row["i"]), int(row["j"])
        value, oxygen, clash = (float(row[c]) for c in ["energy", "oxygen", "clash"])
        if kind == "base_base":
            add(i, j, "E_bb", value)
        elif kind == "base_clash":
            add(i, j, "E_clash", value)
        elif kind == "ribose_clash":
            add(i, j, "E_oo", oxygen)
            add(i, j, "E_clash", clash)
        elif kind == "phosphate_clash":
            add(owner(i), owner(j), "E_clash", value)
        elif kind in {"base_ribose", "base_phosphate", "ribose_phosphate"}:
            target = j if kind == "base_ribose" else owner(j)
            local_column = "E_" + kind
            add(
                i,
                target,
                "E_oo" if kind == "ribose_phosphate" else "E_bo",
                oxygen,
                local_column,
            )
            add(i, target, "E_clash", clash, local_column)
        elif kind == "rotamer":
            local[i]["E_rotamer"] += value
        elif kind == "phosphate_local":
            local[i]["E_phosphate_local"] += value
    for parts in pairs.values():
        parts["E_pair"] = math.fsum(parts[c] for c in PAIR_NUMERIC[:-1])
    for parts in local:
        parts["E_nonpair"] = math.fsum(parts[c] for c in LOCAL_NUMERIC[:-1])
    return pairs, local, phosphate


def verify_export(pair_path, local_path, rows, ids, ownership):
    pairs = read_rows(pair_path)
    local = read_rows(local_path)
    expected_pairs, expected_local, phosphate = aggregate_oracle(
        rows, len(ids), ownership
    )
    assert len(pairs) == len(expected_pairs)
    assert len(local) == len(ids)
    id_to_index = {key: i for i, key in enumerate(ids)}
    assert len(id_to_index) == len(ids), "Duplicate author residue IDs"
    seen = set()
    max_error = 0.0
    for row in pairs:
        i = id_to_index[(row["chain_i"], row["residue_i"])]
        j = id_to_index[(row["chain_j"], row["residue_j"])]
        assert i < j and (i, j) not in seen, "Duplicate/reversed pair"
        seen.add((i, j))
        # Sequence separation is a property of the reported PDB/node endpoints.
        if j == i + 1 and phosphate[i] == i + 1:
            sep = "1"
        elif j == i + 2 and phosphate[i] == i + 1 and phosphate[i + 1] == i + 2:
            sep = "2"
        else:
            sep = "long"
        assert row["sep"] == sep, (i, j, row["sep"], sep)
        for column in PAIR_NUMERIC:
            actual = float(row[column])
            expected = expected_pairs[i, j][column]
            close(actual, expected, f"{ownership} pair {ids[i]} {ids[j]} {column}")
            max_error = max(max_error, abs(actual - expected))
        close(
            float(row["E_pair"]),
            math.fsum(float(row[c]) for c in PAIR_NUMERIC[:-1]),
            "Pair sum",
        )
    assert seen == set(expected_pairs)
    seen_local = set()
    max_nonpair_error = 0.0
    for row in local:
        i = id_to_index[(row["chain"], row["residue"])]
        assert i not in seen_local
        seen_local.add(i)
        for column in LOCAL_NUMERIC:
            close(
                float(row[column]),
                expected_local[i][column],
                f"Local {ids[i]} {column}",
            )
            max_nonpair_error = max(
                max_nonpair_error, abs(float(row[column]) - expected_local[i][column])
            )
    total = math.fsum(float(row["E_pair"]) for row in pairs) + math.fsum(
        float(row["E_nonpair"]) for row in local
    )
    reference = float(next(row["energy"] for row in rows if row["kind"] == "global"))
    close(total, reference, "Pair+nonpair reconstruction")
    return pairs, {
        "pair_count": len(pairs),
        "nonpair_count": len(local),
        "reconstructed": total,
        "reconstruction_error": abs(total - reference),
        "max_pair_component_error": max_error,
        "max_nonpair_component_error": max_nonpair_error,
    }


def compare_raw_components(path, reference):
    raw_rows = read_rows(path, "\t")
    actual = {
        (row["component"], int(row["residue_i"]) - 1, int(row["residue_j"]) - 1): float(
            row["energy"]
        )
        for row in raw_rows
    }
    assert len(actual) == len(raw_rows), "Duplicate directional component"
    expected = {
        (row["kind"], int(row["i"]), int(row["j"])): float(row["energy"])
        for row in reference
        if row["kind"] in RAW_KINDS
    }
    assert set(actual) == set(expected), "Missing/duplicate directional components"
    for key, value in expected.items():
        close(actual[key], value, f"Original helper {key}")
    return {
        "record_count": len(actual),
        "max_error": max(abs(actual[key] - value) for key, value in expected.items()),
    }


def representatives(pairs, pdb):
    selected = {}
    coordinates = defaultdict(list)
    for line in pdb.read_text().splitlines():
        if (
            line.startswith(("ATOM  ", "HETATM"))
            and len(line) >= 54
            and line[16] in " A1"
        ):
            atom = line[12:16].strip()
            if atom and atom[0] not in "HD" and line[17:20].strip() != "HOH":
                coordinates[(line[21:22], line[22:27].strip())].append(
                    tuple(float(line[a:b]) for a, b in [(30, 38), (38, 46), (46, 54)])
                )

    def distance(row):
        a = coordinates[row["chain_i"], row["residue_i"]]
        b = coordinates[row["chain_j"], row["residue_j"]]
        return min(math.dist(x, y) for x in a for y in b)

    for label in ["stacking", "base_pair", "non_watson_crick_pair"]:
        row = next((r for r in pairs if label in r["contact_type"].split(";")), None)
        if row:
            selected[label] = row
    row = next((r for r in pairs if float(r["E_clash"]) > 1e-6), None)
    if row:
        selected["positive_clash"] = row
    row = next((r for r in pairs if r["sep"] == "long" and distance(r) > 12), None)
    if row:
        selected["spatially_separated"] = dict(
            row, minimum_heavy_atom_distance=distance(row)
        )
        close(float(row["E_pair"]), 0, "Spatially separated control")
    return selected


def check_case(name, fixtures, work, runner, image, reference):
    case = work / name
    case.mkdir()
    pdb = case / f"{name}.pdb"
    shutil.copyfile(fixtures / f"{name}.pdb", pdb)
    ids = pdb_ids(pdb)
    template = (fixtures / f"{name}.briq").read_text().splitlines()
    config = case / "input.briq"
    config.write_text(
        "\n".join(
            ["pdb /work/" + pdb.relative_to(work).as_posix()]
            + [line for line in template if not line.startswith("pdb ")]
        )
        + "\n",
        newline="\n",
    )
    prepared = case / "prepared.briq"
    runner.run(image, "BRiQ_Energy", "--prepare", pdb, prepared)
    prepared_lines = prepared.read_text().splitlines()
    assert next(x for x in prepared_lines if x.startswith("seq ")) == next(
        x for x in template if x.startswith("seq ")
    )
    if name == "1RNA":
        assert "break 13" in prepared_lines
        check_chain_boundary(pdb, runner, image)
    raw = case / "original.tsv"
    runner.run(reference, "BRiQ_Oracle", config, raw)
    oracle = read_rows(raw, "\t")
    original_energy = float(next(r["energy"] for r in oracle if r["kind"] == "global"))
    close(original_energy, EXPECTED[name], f"{name} original baseline")
    pair_path = case / "pairs.csv"
    local_path = case / "nonpair.csv"
    components = case / "components.tsv"
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
    )
    close(last_energy(result), original_energy, "Export vs original global")
    pairs, report = verify_export(pair_path, local_path, oracle, ids, "pdb")
    check_local_baselines(name, pairs)
    report["fixed_local_baselines"] = len(LOCAL_BASELINES[name])
    report["directional_components"] = compare_raw_components(components, oracle)
    report["original_global"] = original_energy
    report["exported_global"] = last_energy(result)
    report["global_error"] = abs(last_energy(result) - original_energy)
    report["representatives"] = representatives(pairs, pdb)
    assert {"stacking", "base_pair", "positive_clash", "spatially_separated"} <= set(
        report["representatives"]
    ), "Missing representative interaction category"
    report["source_sha256"] = hashlib.sha256(pdb.read_bytes()).hexdigest()
    report["oxygen_query_reasons"] = dict(
        Counter(r["reason"] for r in read_rows(Path(str(raw) + ".oxygen.tsv"), "\t"))
    )
    report["hbond_query_reasons"] = {
        kind: dict(
            Counter(
                r["reason"]
                for r in read_rows(Path(str(raw) + ".hbond.tsv"), "\t")
                if r["kind"] == kind
            )
        )
        for kind in ["O2O2", "O2OP"]
    }
    report["parameters"] = json.loads(Path(str(raw) + ".parameters.json").read_text())
    report["nonzero_columns"] = {
        c: sum(float(r[c]) != 0 for r in pairs) for c in PAIR_NUMERIC[:-1]
    }
    if name == "1RNA":
        assert report["pair_count"] == 378 and report["nonpair_count"] == 28
        assert (
            report["nonzero_columns"]["E_bo"] == 0
            and report["nonzero_columns"]["E_oo"] == 0
        )
    repeat_pairs = case / "repeat.csv"
    repeat = runner.run(image, "BRiQ_Energy", config, "--pair-energies", repeat_pairs)
    close(last_energy(repeat), original_energy, "Repeated score")
    assert pair_path.read_bytes() == repeat_pairs.read_bytes(), (
        "Pair export is nondeterministic"
    )
    assert (
        local_path.read_bytes() == Path(str(repeat_pairs) + ".nonpair.csv").read_bytes()
    ), "Local export is nondeterministic"
    legacy = runner.run(image, "BRiQ_Energy", config)
    close(last_energy(legacy), original_energy, "Legacy one-argument scoring")
    node_path = case / "node.csv"
    node_result = runner.run(
        image,
        "BRiQ_Energy",
        config,
        "--pair-energies",
        node_path,
        "--pair-ownership",
        "node",
    )
    close(last_energy(node_result), original_energy, "Legacy node ownership score")
    node_pairs, node_report = verify_export(
        node_path, Path(str(node_path) + ".nonpair.csv"), oracle, ids, "node"
    )
    report["node_ownership"] = node_report
    report["reassigned_pair_count"] = sum(
        any(abs(float(a[c]) - float(b[c])) > TOL for c in PAIR_NUMERIC)
        for a, b in zip(pairs, node_pairs)
    )
    weighted_raw = case / "weighted-original.tsv"
    weights = ["0.5", "1.25", "0.75"]
    runner.run(reference, "BRiQ_Oracle", config, weighted_raw, *weights)
    weighted_pairs = case / "weighted.csv"
    weighted_result = runner.run(
        image,
        "BRiQ_Energy",
        config,
        "--pair-energies",
        weighted_pairs,
        "--weights",
        *weights,
        "--seed",
        "0",
    )
    weighted_reference = read_rows(weighted_raw, "\t")
    close(
        last_energy(weighted_result),
        float(next(r["energy"] for r in weighted_reference if r["kind"] == "global")),
        "Weighted original score",
    )
    _, report["weighted"] = verify_export(
        weighted_pairs,
        Path(str(weighted_pairs) + ".nonpair.csv"),
        weighted_reference,
        ids,
        "pdb",
    )
    report["weighted"]["seed"] = 0
    repeated_weighted = case / "weighted-repeat.csv"
    weighted_repeat = runner.run(
        image,
        "BRiQ_Energy",
        config,
        "--pair-energies",
        repeated_weighted,
        "--weights",
        *weights,
        "--seed",
        "0",
    )
    close(
        last_energy(weighted_repeat),
        last_energy(weighted_result),
        "Seeded weighted score changed",
    )
    assert weighted_pairs.read_bytes() == repeated_weighted.read_bytes(), (
        "Seeded pair export is nondeterministic"
    )
    assert (
        Path(str(weighted_pairs) + ".nonpair.csv").read_bytes()
        == Path(str(repeated_weighted) + ".nonpair.csv").read_bytes()
    ), "Seeded nonpair export is nondeterministic"
    if name == "1RNA":
        # Repeat author numbers across chains and include an insertion code.
        changed = []
        remap = {
            key: (key[0], str(i % 14 + 1) + ("A" if i == 15 else ""))
            for i, key in enumerate(ids)
        }
        for line in pdb.read_text().splitlines():
            key = (line[21:22], line[22:27].strip())
            if line.startswith(("ATOM  ", "HETATM")) and key in remap:
                new = remap[key][1]
                number = new.rstrip("A")
                insertion = "A" if new.endswith("A") else " "
                line = line[:22] + f"{int(number):4d}" + insertion + line[27:]
            changed.append(line)
        renumbered = case / "renumbered.pdb"
        renumbered.write_text("\n".join(changed) + "\n", newline="\n")
        alternate = case / "renumbered.briq"
        alternate.write_text(
            config.read_text().replace(f"{name}/{name}.pdb", f"{name}/renumbered.pdb"),
            newline="\n",
        )
        renum_pairs = case / "renumbered.csv"
        renum = runner.run(
            image, "BRiQ_Energy", alternate, "--pair-energies", renum_pairs
        )
        close(
            last_energy(renum),
            original_energy,
            "Renumbering/insertion code changed score",
        )
        verify_export(
            renum_pairs,
            Path(str(renum_pairs) + ".nonpair.csv"),
            oracle,
            pdb_ids(renumbered),
            "pdb",
        )
        report["renumbering_and_insertion_codes"] = "passed"
    print(
        f"PASS {name}: E={original_energy:.17g}, pairs={len(pairs)}, max local error={report['max_pair_component_error']:.3g}",
        flush=True,
    )
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--docker", required=True, help="Candidate image")
    parser.add_argument("--reference-image", default="briq:validation-original")
    parser.add_argument("--build-reference", action="store_true")
    parser.add_argument(
        "--fixtures", type=Path, default=Path(__file__).resolve().parent / "fixtures"
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    work = (
        args.output.resolve()
        if args.output
        else Path(tempfile.mkdtemp(prefix="briq-final-validation-"))
    )
    work.mkdir(parents=True, exist_ok=True)
    assert not (work / "1RNA").exists(), (
        "Use a new output directory to preserve previous results"
    )
    report = {
        "reference_commit": REFERENCE_COMMIT,
        "tolerance": TOL,
        "candidate_image": args.docker,
        "cases": {},
    }
    print("Validation artifacts:", work, flush=True)
    if args.build_reference:
        build_reference(
            Path(__file__).resolve().parents[1], work, args.docker, args.reference_image
        )
    report["reference_image"] = args.reference_image
    report["image_ids"] = {
        image: subprocess.check_output(
            ["docker", "image", "inspect", image, "--format", "{{.Id}}"], text=True
        ).strip()
        for image in [args.docker, args.reference_image]
    }
    runner = Runner(work)
    for name in EXPECTED:
        report["cases"][name] = check_case(
            name,
            args.fixtures.resolve(),
            work,
            runner,
            args.docker,
            args.reference_image,
        )
        (work / "report.json").write_text(json.dumps(report, indent=2))
    assert "non_watson_crick_pair" in report["cases"]["1EHZ"]["representatives"]
    assert report["cases"]["1EHZ"]["nonzero_columns"]["E_bo"] > 0, (
        "No positive control for base-oxygen export"
    )
    assert report["cases"]["1EHZ"]["nonzero_columns"]["E_oo"] > 0, (
        "No positive control for oxygen-oxygen export"
    )
    print(
        "All independent global, directional, pair, ownership and regression checks passed.",
        flush=True,
    )


if __name__ == "__main__":
    main()
