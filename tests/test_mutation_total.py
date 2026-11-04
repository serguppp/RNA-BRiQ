"""Compare memory mutations to independent PDB-batch results and test rejection."""

import argparse
import csv
import json
import math
import shutil
import subprocess
import tempfile
from itertools import product
from pathlib import Path


def read_rows(path):
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def compare(actual, expected):
    keys = ["chain_i", "residue_i", "chain_j", "residue_j", "contact_type", "base_i", "base_j"]
    numeric = ["native_E_total", "E_total", "delta_E_total"]
    assert len(actual) == len(expected), (len(actual), len(expected))
    largest = 0.0
    groups = {}
    for a, b in zip(actual, expected, strict=True):
        assert all(a[k] == b[k] for k in keys), (a, b)
        for key in numeric:
            x, y = float(a[key]), float(b[key])
            assert math.isfinite(x) and math.isfinite(y)
            largest = max(largest, abs(x - y))
            assert abs(x - y) <= 1e-9, (key, a, b)
        pair = tuple(a[k] for k in keys[:4])
        groups.setdefault(pair, []).append((a["base_i"], a["base_j"]))
    for variants in groups.values():
        assert variants == list(product("ACGU", repeat=2))
    return {"contacts": len(groups), "variants": len(actual), "max_difference": largest}


def validate_case(case, expected, image, verify_full=False):
    case = case.resolve()
    summary = compare(read_rows(case / "memory/contact_variants.csv"), read_rows(expected))
    with tempfile.TemporaryDirectory(prefix="briq-memory-validation-") as temporary:
        work = Path(temporary)
        count = 0

        def run(geometry, output, *, ok=False, extra=()):
            nonlocal count
            count += 1
            command = [
                "docker",
                "run",
                "--rm",
                "--mount",
                f"type=bind,source={case},target=/case,readonly",
                "--mount",
                f"type=bind,source={work},target=/work",
                image,
                "BRiQ_Energy",
                "--mutate-total",
                "/case/reference.pdb",
                geometry,
                output,
                *extra,
            ]
            result = subprocess.run(command, capture_output=True, text=True)
            (work / f"test-{count}.log").write_text(result.stdout + result.stderr)
            assert (result.returncode == 0) == ok, (
                command,
                result.stdout[-2000:],
                result.stderr[-2000:],
            )

        original = (case / "geometry.tsv").read_text()
        lines = original.splitlines()
        fields = lines[1].split("\t")
        removed_key = fields[:3]
        corruptions = {
            "header": "invalid\n",
            "empty": lines[0] + "\n",
            "duplicate": original + lines[1] + "\n",
            "nan": original.replace(lines[1], "\t".join(fields[:4] + ["nan", *fields[5:]]), 1),
            "unknown": original.replace(lines[1], "\t".join(["Z", "9999", *fields[2:]]), 1),
            "missing": "\n".join(line for line in lines if line.split("\t")[:3] != removed_key)
            + "\n",
        }
        # Change the first O3' coordinate so the fast path would otherwise retain
        # reference separations while the PDB scorer detects a new break.
        for line in lines[1:]:
            f = line.split("\t")
            if f[3] == "O3'":
                f[4:7] = ["999", "999", "999"]
                corruptions["topology"] = original.replace(line, "\t".join(f), 1)
                break
        for name, text in corruptions.items():
            (work / "bad.tsv").write_text(text, newline="\n")
            run("/work/bad.tsv", f"/work/rejected-{name}")
            assert not (work / f"rejected-{name}").exists(), name
        occupied = work / "occupied"
        occupied.mkdir()
        sentinel = occupied / "keep.txt"
        sentinel.write_text("preserve existing results")
        run("/case/geometry.tsv", "/work/occupied")
        assert sentinel.read_text() == "preserve existing results"
        assert list(occupied.iterdir()) == [sentinel]
        if verify_full:
            run("/case/geometry.tsv", "/work/verified", ok=True, extra=("--verify-full",))
            compare(read_rows(work / "verified/contact_variants.csv"), read_rows(expected))
            assert json.loads((work / "verified/mutation_run.json").read_text())["verify_full"]
        summary["rejection_cases"] = len(corruptions) + 1
        summary["verify_full"] = verify_full
    (case / "validation.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--image", default="briq:memory")
    parser.add_argument("--verify-full", action="store_true")
    args = parser.parse_args()
    fixtures = Path(__file__).resolve().parent / "fixtures/memory"
    with tempfile.TemporaryDirectory(prefix="briq-memory-cases-") as temporary:
        for name in ("1RNA", "1EHZ"):
            work = Path(temporary) / name
            work.mkdir()
            shutil.copyfile(fixtures / f"{name}.pdb", work / "reference.pdb")
            shutil.copyfile(fixtures / f"{name}.tsv", work / "geometry.tsv")
            with (work / "memory.log").open("w") as log:
                subprocess.run(
                    [
                        "docker",
                        "run",
                        "--rm",
                        "--mount",
                        f"type=bind,source={work},target=/work",
                        args.image,
                        "BRiQ_Energy",
                        "--mutate-total",
                        "/work/reference.pdb",
                        "/work/geometry.tsv",
                        "/work/memory",
                    ],
                    check=True,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                )
            validate_case(work, fixtures / f"{name}_expected.csv", args.image, args.verify_full)


if __name__ == "__main__":
    main()
