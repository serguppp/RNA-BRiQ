"""Verify checkpoint resume, legacy-prefix adoption and rejection before append."""

import argparse
import csv
import shutil
import subprocess
import tempfile
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--image", default="briq:batch")
    args = parser.parse_args()
    fixtures = Path(__file__).resolve().parent / "fixtures"
    with tempfile.TemporaryDirectory(prefix="briq-resume-") as directory:
        work = Path(directory)
        shutil.copyfile(fixtures / "1RNA.pdb", work / "native, input.pdb")
        shutil.copyfile(fixtures / "1EHZ.pdb", work / "1EHZ.pdb")
        # Repeated input after a different length also tests cache reset.
        paths = ["/work/native, input.pdb", "/work/1EHZ.pdb", "/work/native, input.pdb"]
        manifest = work / "jobs.tsv"
        manifest.write_text("\n".join(paths) + "\n")
        counter = 0

        def run(output, *options, ok=True):
            nonlocal counter
            counter += 1
            log = work / f"log-{counter}.txt"
            with log.open("w") as stream:
                result = subprocess.run(
                    [
                        "docker",
                        "run",
                        "--rm",
                        "--mount",
                        f"type=bind,source={work},target=/work",
                        args.image,
                        "BRiQ_Energy",
                        "--batch-total",
                        "/work/jobs.tsv",
                        f"/work/{output}",
                        "--engine",
                        "incremental",
                        *options,
                    ],
                    stdout=stream,
                    stderr=subprocess.STDOUT,
                    check=False,
                )
            text = log.read_text(errors="replace")
            assert (result.returncode == 0) == ok, text[-12000:]
            return text

        run("fresh.csv", "--flush-every", "1")
        complete = (work / "fresh.csv").read_bytes()
        lines = complete.splitlines(keepends=True)
        state = (work / "fresh.csv.state").read_bytes()
        assert len(lines) == 4
        with (work / "fresh.csv").open(newline="") as file:
            rows = list(csv.DictReader(file))
        assert rows[0]["E_total"] == rows[2]["E_total"]
        print("Fresh batch and repeated-input cache reset PASS", flush=True)

        (work / "fresh.csv").write_bytes(b"".join(lines[:3]))
        trusted_state = b"\n".join(state.splitlines()[:-1]) + b"\n"
        (work / "fresh.csv.state").write_bytes(trusted_state)
        log = run("fresh.csv", "--resume")
        assert (
            "[batch-total 3/3]" in log
            and "[batch-total 1/3]" not in log
            and "[batch-total 2/3]" not in log
        )
        assert (work / "fresh.csv").read_bytes() == complete
        print("Committed prefix skipped without recomputing PASS", flush=True)

        (work / "fresh.csv").write_bytes(b"".join(lines[:3]))
        interrupted_state = b"\n".join(state.splitlines()[:-2]) + b"\n"
        (work / "fresh.csv.state").write_bytes(interrupted_state)
        log = run("fresh.csv", "--resume")
        assert "[batch-total 2/3]" in log and "[batch-total 1/3]" not in log
        assert (work / "fresh.csv").read_bytes() == complete
        print(
            "Checkpoint resume preserves prefix and verifies uncommitted tail PASS",
            flush=True,
        )

        log = run("fresh.csv", "--resume")
        assert "Loading shared energy tables" not in log
        assert (work / "fresh.csv").read_bytes() == complete

        def rejected(contents, checkpoint=state):
            target = work / "reject.csv"
            target.write_bytes(contents)
            (work / "reject.csv.state").write_bytes(checkpoint)
            log = run("reject.csv", "--resume", ok=False)
            assert "Loading shared energy tables" not in log
            assert target.read_bytes() == contents

        rejected(complete.replace(lines[1], lines[1].replace(b"-28.", b"99928.")))
        rejected(b"bad header\n")
        rejected(lines[0] + lines[2])
        rejected(lines[0] + lines[1].replace(b"-28.", b"nan28."))
        rejected(b"".join(lines[:2]).rstrip(b"\n"))
        rejected(b"".join(lines[:2]), state + b"tampered\n")
        original = (work / "native, input.pdb").read_bytes()
        (work / "native, input.pdb").write_bytes(original + b"REMARK changed\n")
        rejected(b"".join(lines[:2]))
        (work / "native, input.pdb").write_bytes(original)
        run("fresh.csv", ok=False)
        print(
            "Input/source identity, invalid CSV and overwrite rejection PASS",
            flush=True,
        )

        legacy = work / "legacy.csv"
        legacy.write_bytes(b"".join(lines[:3]))
        run("legacy.csv", "--resume", ok=False)
        assert legacy.read_bytes() == b"".join(lines[:3])
        run("legacy.csv", "--resume-verify")
        assert legacy.read_bytes() == complete
        assert (work / "legacy.csv.state").is_file()
        print(
            "Legacy CSV adoption verifies all saved energies before append PASS",
            flush=True,
        )

        wrong = work / "wrong.csv"
        wrong_lines = list(lines[:2])
        fields = next(csv.reader([wrong_lines[1].decode().strip()]))
        wrong_lines[1] = wrong_lines[1].replace(fields[2].encode(), b"999")
        wrong.write_bytes(b"".join(wrong_lines))
        run("wrong.csv", "--resume-verify", ok=False)
        assert wrong.read_bytes() == b"".join(wrong_lines)
        assert not (work / "wrong.csv.state").exists()
        print(
            "Mismatched legacy energy rejected without modifying saved rows PASS",
            flush=True,
        )


if __name__ == "__main__":
    main()
