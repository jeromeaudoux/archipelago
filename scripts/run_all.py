"""Run the full precompute pipeline: build_bundles -> add_clinvar -> build_index.

Usage:
  python3 scripts/run_all.py                 # all MANE genes
  python3 scripts/run_all.py --genes KIF1A,TP53
"""
import argparse
import subprocess
import sys
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def run(script, extra):
    cmd = [sys.executable, os.path.join(HERE, script)] + extra
    print(f"\n=== {script} {' '.join(extra)} ===", flush=True)
    subprocess.run(cmd, check=True, cwd=HERE)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes")
    args = ap.parse_args()
    import os, config as C
    extra = ["--genes", args.genes] if args.genes else []
    if not os.path.exists(C.UNIPROT_MAP_JSON):
        run("build_uniprot_map.py", [])
    run("build_bundles.py", extra)
    run("add_clinvar.py", extra)
    run("build_exons.py", [])
    run("add_hotspots.py", extra)
    run("add_rmc.py", extra)
    run("build_index.py", [])


if __name__ == "__main__":
    main()
