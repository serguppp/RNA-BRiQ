"""Integration test for the five-file BRiQ output mode (Docker + stdlib)."""

import argparse
import csv
import json
import math
import subprocess
import tempfile
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--docker", default="briq:output-dir")
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    fixtures = root / "fixtures"
    with tempfile.TemporaryDirectory(prefix="briq-output-") as tmp:
        work = Path(tmp)
        base = [
            "docker",
            "run",
            "--rm",
            "--mount",
            f"type=bind,source={fixtures},target=/fixtures,readonly",
            "--mount",
            f"type=bind,source={work},target=/work",
            args.docker,
            "BRiQ_Energy",
        ]

        def run(*options, ok=True):
            result = subprocess.run(
                base + list(options), check=False, capture_output=True, text=True
            )
            assert (result.returncode == 0) == ok, result.stdout + result.stderr
            return result

        for name, expected in [
            ("1RNA", -28.606109849613684),
            ("1EHZ", -133.66638013285746),
        ]:
            expected_files = {
                f"{name}.pdb",
                f"{name}.briq",
                "pairs.csv",
                "nonpair.csv",
                "run.json",
            }
            out = work / name
            run(f"/fixtures/{name}.pdb", "--output-dir", f"/work/{name}", "--seed", "0")
            assert {p.name for p in out.iterdir()} == expected_files
            assert (out / f"{name}.pdb").read_bytes() == (
                fixtures / f"{name}.pdb"
            ).read_bytes()
            assert (out / f"{name}.briq").read_text().startswith(f"pdb {name}.pdb\n")
            meta = json.loads((out / "run.json").read_text())
            assert meta["program"] == "BRiQ_Energy"
            assert meta["source_version"] != "unknown"
            assert meta["parameters"] == {
                "break_weight": 1,
                "connect_weight": 1,
                "clash_weight": 1,
                "seed": 0,
                "pair_ownership": "pdb",
            }
            assert math.isclose(meta["energy_total"], expected, abs_tol=1e-9)
            total = 0
            for file, column in [("pairs.csv", "E_pair"), ("nonpair.csv", "E_nonpair")]:
                with (out / file).open(newline="") as handle:
                    total += math.fsum(
                        float(row[column]) for row in csv.DictReader(handle)
                    )
            assert math.isclose(total, meta["energy_total"], abs_tol=1e-9)
            before = {p.name: p.read_bytes() for p in out.iterdir()}
            run(f"/fixtures/{name}.pdb", "--output-dir", f"/work/{name}", ok=False)
            assert {p.name: p.read_bytes() for p in out.iterdir()} == before
            print(name, "PASS", meta["energy_total"], flush=True)

        run(
            "/fixtures/1RNA.pdb",
            "--output-dir",
            "/work/conflict",
            "--components",
            "/work/components.tsv",
            ok=False,
        )
        assert not (work / "conflict").exists()
        run("/fixtures/missing.pdb", "--output-dir", "/work/missing", ok=False)
        assert not (work / "missing").exists()
        run("/fixtures/1RNA.pdb", "--output-dir", ok=False)
        print("Argument validation and overwrite protection PASS", flush=True)


if __name__ == "__main__":
    main()
