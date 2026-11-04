"""Validate total-only batches against full exports (Docker + stdlib)."""

import argparse
import csv
import json
import math
import shutil
import subprocess
import tempfile
import time
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--docker", default="briq:batch")
    parser.add_argument(
        "--variants", type=Path, help="Optional saved 1RNA variant scores"
    )
    args = parser.parse_args()
    fixtures = Path(__file__).resolve().parent / "fixtures"
    with tempfile.TemporaryDirectory(prefix="briq-total-") as directory:
        work = Path(directory)
        scratch = work / "scratch"
        scratch.mkdir()
        # The original PDB names need not be globally unique.
        special = work / "input with spaces,comma.pdb"
        shutil.copyfile(fixtures / "1RNA.pdb", special)
        inputs = [
            "/fixtures/1RNA.pdb",
            "/fixtures/1EHZ.pdb",
            "/fixtures/1RNA.pdb",
            "/work/input with spaces,comma.pdb",
        ]
        if args.variants:
            for base in "ACGU":
                name = f"1RNA_A1{base}"
                source = args.variants / name
                shutil.copyfile(source / f"{name}.pdb", work / f"{name}.pdb")
                inputs.append(f"/work/{name}.pdb")
        mounts = [
            "--mount",
            f"type=bind,source={fixtures},target=/fixtures,readonly",
            "--mount",
            f"type=bind,source={work},target=/work",
            "--mount",
            f"type=bind,source={scratch},target=/tmp",
        ]
        counter = 0

        def run(*options, ok=True):
            nonlocal counter
            counter += 1
            log = work / f"command-{counter}.log"
            start = time.perf_counter()
            with log.open("w") as stream:
                result = subprocess.run(
                    [
                        "docker",
                        "run",
                        "--rm",
                        *mounts,
                        args.docker,
                        "BRiQ_Energy",
                        *options,
                    ],
                    stdout=stream,
                    stderr=subprocess.STDOUT,
                    check=False,
                )
            output = log.read_text(errors="replace")
            assert (result.returncode == 0) == ok, output[-12000:]
            assert not list(scratch.iterdir()), "Temporary preparation files leaked"
            return time.perf_counter() - start, output

        def manifest(text):
            (work / "total.tsv").write_text(text, encoding="utf-8", newline="")

        for text in [
            "",
            "1RNA.pdb\n",
            "/fixtures/1RNA.pdb\t/work/out\n",
            "/fixtures/missing.pdb\n",
            "/fixtures/1RNA.pdb\n/fixtures/missing.pdb\n",
        ]:
            manifest(text)
            _, log = run(
                "--batch-total", "/work/total.tsv", "/work/rejected.csv", ok=False
            )
            assert not (work / "rejected.csv").exists()
            assert "Loading shared energy tables" not in log
        manifest("/fixtures/1RNA.pdb\n")
        marker = work / "keep.csv"
        marker.write_text("keep")
        for output in [
            "/work/keep.csv",
            "/work/missing-parent/out.csv",
            "relative.csv",
            "/work/total.tsv",
            "/fixtures/1RNA.pdb",
        ]:
            _, log = run("--batch-total", "/work/total.tsv", output, ok=False)
            assert "Loading shared energy tables" not in log
        assert marker.read_text() == "keep"
        run("--batch-total", "/work/total.tsv", ok=False)
        print("Preflight, arguments and overwrite protection PASS", flush=True)

        baseline = "".join(f"{path}\t/work/full-{i}\n" for i, path in enumerate(inputs))
        (work / "full.tsv").write_text(baseline)
        full_time, _ = run("--batch", "/work/full.tsv")
        expected = [
            json.loads((work / f"full-{i}" / "run.json").read_text())["energy_total"]
            for i in range(len(inputs))
        ]
        # Historical saved scores predate the phosphate memory-safety fix.
        # The full exporter in the same corrected build is the reference.
        print(
            f"Full-export reference: {len(inputs)} jobs, {full_time:.3f}s", flush=True
        )

        before = {p.relative_to(work) for p in work.rglob("*")}
        # Accept Windows manifests and blank rows; preserve repeated inputs and order.
        manifest("\r\n" + "\r\n".join(inputs) + "\r\n")
        total_time, log = run("--batch-total", "/work/total.tsv", "/work/totals.csv")
        assert log.count("Loading shared energy tables") == 1
        assert "source_version=unknown" not in log
        with (work / "totals.csv").open(newline="") as stream:
            reader = csv.DictReader(stream)
            assert reader.fieldnames == ["job_index", "input_pdb", "E_total"]
            rows = list(reader)
        assert len(rows) == len(inputs)
        for i, (row, path, reference) in enumerate(
            zip(rows, inputs, expected, strict=True)
        ):
            assert int(row["job_index"]) == i
            assert row["input_pdb"] == path
            assert math.isfinite(float(row["E_total"]))
            assert math.isclose(
                float(row["E_total"]), reference, rel_tol=0, abs_tol=1e-9
            )
        assert rows[0]["E_total"] == rows[2]["E_total"] == rows[3]["E_total"]
        after = {p.relative_to(work) for p in work.rglob("*")}
        created = after - before
        assert created == {
            Path("totals.csv"),
            Path("totals.csv.state"),
            Path(f"command-{counter}.log"),
        }
        print(
            "Energies, order, CSV escaping, repeatability and temporary cleanup PASS",
            flush=True,
        )
        print(
            json.dumps(
                {
                    "jobs": len(inputs),
                    "full_batch_seconds": full_time,
                    "total_batch_seconds": total_time,
                    "speedup": full_time / total_time,
                }
            ),
            flush=True,
        )


if __name__ == "__main__":
    main()
