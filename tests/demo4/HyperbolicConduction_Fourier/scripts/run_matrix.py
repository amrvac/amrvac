#!/usr/bin/env python3
"""Run the generated matrix sequentially, stopping on the first failure."""

from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TIMELOOP_RE = re.compile(r"Total timeloop took\s*:\s*([0-9.Ee+-]+) sec")
FINISHED_RE = re.compile(r"Finished AMRVAC in\s*:\s*([0-9.Ee+-]+) sec")


def load_rows(path: Path, suite: str, cases: set[str]) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if suite != "all":
        rows = [row for row in rows if row["suite"] == suite]
    if cases:
        rows = [row for row in rows if row["case"] in cases]
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=Path, default=ROOT / "matrix.csv")
    parser.add_argument("--suite", choices=("all", "aniso", "isotropic"), default="all")
    parser.add_argument("--case", action="append", default=[])
    parser.add_argument("--ranks", type=int, default=4)
    parser.add_argument("--timeout", type=float, default=1800.0)
    parser.add_argument("--resume", action="store_true", help="Skip only cases with successful metadata and outputs")
    args = parser.parse_args()

    rows = load_rows(args.matrix, args.suite, set(args.case))
    if not rows:
        raise SystemExit("No matrix rows selected")
    executable = ROOT / "amrvac"
    if not executable.is_file():
        raise SystemExit(f"Missing executable: {executable}")

    print(f"Selected {len(rows)} cases; no failed case will be retried automatically", flush=True)
    for index, row in enumerate(rows, start=1):
        run_dir = ROOT / row["run_dir"]
        metadata_path = run_dir / "run_metadata.json"
        snapshot = run_dir / "output0000.dat"
        log_path = run_dir / "output.log"
        if args.resume and metadata_path.exists() and snapshot.exists() and log_path.exists():
            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
            if metadata.get("returncode") == 0:
                print(f"[{index:02d}/{len(rows):02d}] SKIP verified {row['case']}", flush=True)
                continue
        if any(run_dir.glob("output*")) or metadata_path.exists():
            raise SystemExit(
                f"Refusing to overwrite existing incomplete/unverified run {run_dir}; inspect it explicitly"
            )
        run_dir.mkdir(parents=True, exist_ok=True)
        command = [
            "mpirun",
            "-np",
            str(args.ranks),
            str(executable),
            "-i",
            "common.par",
            row["parfile"],
            row["overlay"],
        ]
        print(f"[{index:02d}/{len(rows):02d}] START {row['case']}", flush=True)
        started = datetime.now(timezone.utc).isoformat()
        tic = time.perf_counter()
        try:
            result = subprocess.run(
                command,
                cwd=ROOT,
                text=True,
                capture_output=True,
                timeout=args.timeout,
                check=False,
            )
            timed_out = False
        except subprocess.TimeoutExpired as exc:
            wall = time.perf_counter() - tic
            (run_dir / "stdout.txt").write_text(exc.stdout or "", encoding="utf-8")
            (run_dir / "stderr.txt").write_text(exc.stderr or "", encoding="utf-8")
            metadata_path.write_text(
                json.dumps(
                    {
                        **row,
                        "command": command,
                        "started_utc": started,
                        "wall_clock_s": wall,
                        "returncode": None,
                        "timed_out": True,
                    },
                    indent=2,
                ),
                encoding="utf-8",
            )
            raise SystemExit(f"TIMEOUT {row['case']} after {wall:.1f} s; matrix stopped without retry")

        wall = time.perf_counter() - tic
        stdout = result.stdout or ""
        stderr = result.stderr or ""
        (run_dir / "stdout.txt").write_text(stdout, encoding="utf-8")
        (run_dir / "stderr.txt").write_text(stderr, encoding="utf-8")
        timeloop = TIMELOOP_RE.search(stdout)
        finished = FINISHED_RE.search(stdout)
        metadata = {
            **row,
            "command": command,
            "started_utc": started,
            "completed_utc": datetime.now(timezone.utc).isoformat(),
            "wall_clock_s": wall,
            "amrvac_timeloop_s": float(timeloop.group(1)) if timeloop else None,
            "amrvac_finished_s": float(finished.group(1)) if finished else None,
            "returncode": result.returncode,
            "timed_out": timed_out,
            "snapshot_exists": snapshot.exists(),
            "log_exists": log_path.exists(),
        }
        metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
        if result.returncode != 0 or not snapshot.exists() or not log_path.exists():
            print(stderr[-2000:], flush=True)
            raise SystemExit(
                f"FAILED {row['case']} rc={result.returncode}; matrix stopped without retry"
            )
        print(
            f"[{index:02d}/{len(rows):02d}] DONE  {row['case']} wall={wall:.3f}s "
            f"timeloop={metadata['amrvac_timeloop_s']}",
            flush=True,
        )


if __name__ == "__main__":
    main()
