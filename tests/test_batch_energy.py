"""Compare batch scoring with fresh BRiQ processes (Docker + stdlib)."""

import argparse
import csv
import json
import math
import subprocess
import tempfile
import time
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--docker", default="briq:batch")
    parser.add_argument("--reference", default="briq:output-dir")
    args = parser.parse_args()
    fixtures = Path(__file__).resolve().parent / "fixtures"

    with tempfile.TemporaryDirectory(prefix="briq-batch-") as directory:
        work = Path(directory)
        mounts = [
            "--mount",
            f"type=bind,source={fixtures},target=/fixtures,readonly",
            "--mount",
            f"type=bind,source={work},target=/work",
        ]
        counter = 0

        def run(image, *options, ok=True):
            nonlocal counter
            counter += 1
            start = time.perf_counter()
            with (work / f"command-{counter}.log").open("w") as log:
                result = subprocess.run(
                    ["docker", "run", "--rm", *mounts, image, "BRiQ_Energy", *options],
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    check=False,
                )
            elapsed = time.perf_counter() - start
            output = (work / f"command-{counter}.log").read_text(errors="replace")
            assert (result.returncode == 0) == ok, output[-12000:]
            return elapsed, output

        def manifest(text):
            (work / "jobs.tsv").write_text(text, encoding="utf-8", newline="")

        invalid = [
            ("", "empty"),
            ("/fixtures/1RNA.pdb /work/no-output\n", "expected PDB"),
            ("1RNA.pdb\t/work/no-output\n", "absolute"),
            ("/fixtures/missing.pdb\t/work/no-output\n", "Cannot read"),
            ("/fixtures/1RNA.pdb\t/work/missing-parent/out\n", "parent"),
            (
                "/fixtures/1RNA.pdb\t/work/no-output\n"
                "/fixtures/missing.pdb\t/work/later\n",
                "Cannot read",
            ),
            (
                "/fixtures/1RNA.pdb\t/work/no-output\n"
                "/fixtures/1EHZ.pdb\t/work/no-output\n",
                "Duplicate",
            ),
        ]
        for text, message in invalid:
            manifest(text)
            _, output = run(args.docker, "--batch", "/work/jobs.tsv", ok=False)
            assert message in output, output
            assert not (work / "no-output").exists()
            assert "Loading shared energy tables" not in output
        occupied = work / "occupied"
        occupied.mkdir()
        (occupied / "keep.txt").write_text("unchanged")
        manifest("/fixtures/1RNA.pdb\t/work/occupied\n")
        _, output = run(args.docker, "--batch", "/work/jobs.tsv", ok=False)
        assert "must be empty" in output
        assert (occupied / "keep.txt").read_text() == "unchanged"
        assert list(occupied.iterdir()) == [occupied / "keep.txt"]
        print(
            "Manifest validation, preflight and overwrite protection PASS", flush=True
        )

        names = ["1RNA", "1EHZ", "1RNA"]
        single_seconds = 0.0
        for index, name in enumerate(names):
            elapsed, _ = run(
                args.reference,
                f"/fixtures/{name}.pdb",
                "--output-dir",
                f"/work/single-{index}",
                "--seed",
                "0",
            )
            single_seconds += elapsed
            print(
                f"Reference {index + 1}/{len(names)}: {name} ({elapsed:.3f}s)",
                flush=True,
            )

        # CRLF, a blank row, spaces, an existing empty directory and repeated input.
        (work / "batch 0").mkdir()
        manifest(
            "\r\n"
            + "".join(
                f"/fixtures/{name}.pdb\t/work/batch {index}\r\n"
                for index, name in enumerate(names)
            )
        )
        batch_seconds, output = run(args.docker, "--batch", "/work/jobs.tsv")
        assert output.count("Loading shared energy tables") == 1, output
        assert "[batch 3/3]" in output

        for index, name in enumerate(names):
            reference = work / f"single-{index}"
            candidate = work / f"batch {index}"
            expected = {
                f"{name}.pdb",
                f"{name}.briq",
                "pairs.csv",
                "nonpair.csv",
                "run.json",
            }
            assert {p.name for p in candidate.iterdir()} == expected
            assert (candidate / f"{name}.pdb").read_bytes() == (
                fixtures / f"{name}.pdb"
            ).read_bytes()
            assert (candidate / f"{name}.briq").read_bytes() == (
                reference / f"{name}.briq"
            ).read_bytes()
            meta = json.loads((candidate / "run.json").read_text())
            original_meta = json.loads((reference / "run.json").read_text())
            assert meta["source_version"] != "unknown"
            assert meta["exporter_version"] == "0.3.0"
            assert meta["parameters"] == original_meta["parameters"]
            assert math.isfinite(meta["energy_total"])
            assert math.isclose(
                meta["energy_total"],
                original_meta["energy_total"],
                rel_tol=0,
                abs_tol=1e-9,
            )
            reconstructed = 0.0
            for filename, total_column in [
                ("pairs.csv", "E_pair"),
                ("nonpair.csv", "E_nonpair"),
            ]:
                with (reference / filename).open(newline="") as file:
                    expected_rows = list(csv.DictReader(file))
                with (candidate / filename).open(newline="") as file:
                    actual_rows = list(csv.DictReader(file))
                assert len(expected_rows) == len(actual_rows)
                for expected_row, actual_row in zip(
                    expected_rows, actual_rows, strict=True
                ):
                    assert expected_row.keys() == actual_row.keys()
                    for key in expected_row:
                        if key.startswith("E_"):
                            assert math.isclose(
                                float(expected_row[key]),
                                float(actual_row[key]),
                                rel_tol=0,
                                abs_tol=1e-9,
                            ), (name, filename, key, expected_row, actual_row)
                        else:
                            assert expected_row[key] == actual_row[key]
                reconstructed += math.fsum(
                    float(row[total_column]) for row in actual_rows
                )
            assert math.isclose(
                reconstructed, meta["energy_total"], rel_tol=0, abs_tol=1e-9
            )
            print(name, "all energy components PASS", meta["energy_total"], flush=True)

        for filename in ["pairs.csv", "nonpair.csv"]:
            assert (work / "batch 0" / filename).read_bytes() == (
                work / "batch 2" / filename
            ).read_bytes()
        print("Repeated input after a different structure PASS", flush=True)
        print(
            json.dumps(
                {
                    "jobs": len(names),
                    "single_seconds": single_seconds,
                    "batch_seconds": batch_seconds,
                    "speedup": single_seconds / batch_seconds,
                }
            ),
            flush=True,
        )


if __name__ == "__main__":
    main()
