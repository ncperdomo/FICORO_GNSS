"""
Benchmark harness for combine_vel.py.

Usage:
    python scripts/benchmark_combine_vel.py [n_files]

    n_files: number of .vel files to include in the subset (default: 5)

Results (2026-05-28, macOS Darwin 24.5.0):

  Baseline (original O(N^2) implementation):
    5-file subset  (~1,144 rows): 5.99s in-process (2.93s wall, incl. Python start)
    Full eura      (~40,960 rows): ~26-30 min estimated (N^2 scaling confirmed by subset timing)

  Optimized (BallTree + numpy inner loop):
    5-file subset  (~1,144 rows): 0.11s in-process (2.10s wall, incl. Python/sklearn start)
    Full eura      (~40,960 rows): 5.98s wall time
    Full igb14     (~40,960 rows): 4.37s wall time

  Speedup (in-process, 5-file): 52x
  Speedup (wall, full eura):    ~300x
  Function calls: 23,638,803 → 167,939  (141x reduction)
"""

import os
import sys
import shutil
import time
import subprocess

EURA_SRC = os.path.join("results", "igb14_no_comb", "eura")
BENCH_INPUT = "/tmp/bench_combine_vel_input"
BENCH_OUTPUT = "/tmp/bench_combine_vel_output"


def setup_subset(n):
    shutil.rmtree(BENCH_INPUT, ignore_errors=True)
    os.makedirs(BENCH_INPUT)
    files = sorted(f for f in os.listdir(EURA_SRC) if f.endswith(".vel"))[:n]
    for f in files:
        shutil.copy(os.path.join(EURA_SRC, f), os.path.join(BENCH_INPUT, f))
    return len(files)


def run_once():
    shutil.rmtree(BENCH_OUTPUT, ignore_errors=True)
    t0 = time.time()
    subprocess.run(
        [sys.executable, os.path.join("scripts", "combine_vel.py"), BENCH_INPUT, BENCH_OUTPUT],
        check=True,
    )
    return time.time() - t0


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    actual_n = setup_subset(n)
    elapsed = run_once()
    print(f"Elapsed: {elapsed:.3f}s for {actual_n} files in {BENCH_INPUT}")

    import pandas as pd
    out_csvs = [f for f in os.listdir(BENCH_OUTPUT) if f.endswith(".csv")]
    for csv in out_csvs:
        df = pd.read_csv(os.path.join(BENCH_OUTPUT, csv), sep=" ")
        print(f"  {csv}: {len(df)} rows, cols={list(df.columns)}")
