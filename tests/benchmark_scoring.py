"""Measure scorer engines on one prepared manifest and verify equal energies."""

import argparse
import csv
import json
import math
import statistics
import subprocess
import time
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--manifest", default="jobs-benchmark.tsv")
    parser.add_argument("--image", default="briq:batch")
    parser.add_argument("--engines", nargs="+", default=["legacy", "shared"])
    parser.add_argument("--prefix", default="benchmark")
    args = parser.parse_args()
    work = args.workspace.resolve()
    expected = (work / args.manifest).read_text().splitlines()
    report = {}
    reference = None
    for engine in args.engines:
        output = work / f"{args.prefix}-{engine}.csv"
        log_path = work / f"{args.prefix}-{engine}.log"
        if output.exists():
            raise FileExistsError(output)
        start = time.perf_counter()
        with log_path.open("w") as log:
            subprocess.run(
                [
                    "docker",
                    "run",
                    "--rm",
                    "--mount",
                    f"type=bind,source={work},target=/work",
                    args.image,
                    "BRiQ_Energy",
                    "--batch-total",
                    f"/work/{args.manifest}",
                    f"/work/{output.name}",
                    "--engine",
                    engine,
                    "--profile",
                ],
                check=True,
                stdout=log,
                stderr=subprocess.STDOUT,
            )
        elapsed = time.perf_counter() - start
        with output.open(newline="") as file:
            rows = list(csv.DictReader(file))
        assert len(rows) == len(expected)
        energies = []
        for index, row in enumerate(rows):
            assert int(row["job_index"]) == index
            assert row["input_pdb"] == expected[index]
            energy = float(row["E_total"])
            assert math.isfinite(energy)
            energies.append(energy)
        if reference is None:
            reference = energies
        maximum = max(abs(a - b) for a, b in zip(reference, energies, strict=True))
        assert maximum < 1e-9, (engine, maximum)
        profiles = {}
        resource_seconds = None
        for line in log_path.read_text(errors="replace").splitlines():
            if line.startswith("[profile] "):
                values = [float(v) for v in line.split(" ", 1)[1].split(",")]
                profiles[int(values[0])] = values[1:]
            elif line.startswith("[profile-resources] "):
                resource_seconds = float(line.split("seconds=", 1)[1])
        assert len(profiles) == len(expected)
        groups = {}
        for index, path in enumerate(expected):
            molecule = path.split("/")[2]
            groups.setdefault(molecule, []).append(profiles[index])
        group_stats = {}
        for molecule, samples in groups.items():
            means = [statistics.fmean(row[k] for row in samples) for k in range(5)]
            group_stats[molecule] = dict(
                zip(
                    ["prepare_s", "rotamers_s", "fragments_s", "tree_s", "energy_s"],
                    means,
                    strict=True,
                )
            )
            group_stats[molecule]["mean_job_s"] = sum(means)
        report[engine] = {
            "wall_seconds": elapsed,
            "resources_s": resource_seconds,
            "max_energy_difference": maximum,
            "molecules": group_stats,
        }
        print(engine, json.dumps(report[engine]), flush=True)
    (work / f"{args.prefix}.json").write_text(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
