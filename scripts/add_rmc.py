"""Inject gnomAD regional missense constraint (RMC) into gene bundles.

Reads gnomAD_v2.1.1_transcripts_with_rmc.tsv (per-transcript sub-regions with an
observed/expected missense ratio) and, for genes present, adds an `rmc` array to
the bundle. Regions are placed in protein/residue space; a gene is skipped when
the RMC transcript's length is discordant with the AlphaMissense isoform.

Usage:  python3 add_rmc.py [--genes KIF1A,PTEN]
"""
import argparse
import csv
import gzip
import json
import os
import re
import sys
from collections import defaultdict

import config as C

RESI_RE = re.compile(r"(\d+)$")


def load_rmc():
    """gene -> transcript -> list of regions {s,e,oe,obs,exp,p}."""
    by_gene = defaultdict(lambda: defaultdict(list))
    with open(C.GNOMAD_RMC_TSV) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            gene, tx = r["gene_name"].strip(), r["transcript"].strip()
            ms, me = RESI_RE.search(r["start_aa"] or ""), RESI_RE.search(r["stop_aa"] or "")
            if not (gene and ms and me):
                continue
            try:
                oe = float(r["oe"]); obs = int(float(r["obs"])); exp = float(r["exp"])
                pval = float(r["p_value"])
            except (ValueError, KeyError):
                continue
            a, c = int(ms.group(1)), int(me.group(1))   # genomic order → normalize to residue order
            by_gene[gene][tx].append({
                "s": min(a, c), "e": max(a, c),
                "oe": round(oe, 4), "obs": obs, "exp": round(exp, 1), "p": pval,
            })
    return by_gene


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes")
    args = ap.parse_args()
    targets = set(args.genes.split(",")) if args.genes else None

    if not os.path.isdir(C.GENES_DIR):
        print("No bundles found — run build_bundles.py first.", file=sys.stderr); sys.exit(1)
    print(f"Reading {C.GNOMAD_RMC_TSV}…", flush=True)
    by_gene = load_rmc()
    print(f"  {len(by_gene)} genes with RMC", flush=True)

    updated = skipped_isoform = 0
    for gene, txs in by_gene.items():
        if targets is not None and gene not in targets:
            continue
        path = os.path.join(C.GENES_DIR, gene + ".json.gz")
        if not os.path.exists(path):
            continue
        with gzip.open(path, "rt") as fh:
            b = json.load(fh)
        length = b["length"]
        # Pick the transcript whose region span best fits the AM isoform length.
        best_tx, best_fit = None, None
        for tx, regs in txs.items():
            maxaa = max(r["e"] for r in regs)
            fit = abs(maxaa - length)
            if maxaa <= length + 1 and (best_fit is None or fit < best_fit):
                best_tx, best_fit = tx, fit
        if best_tx is None:
            skipped_isoform += 1
            if "rmc" not in b:
                continue
            del b["rmc"]                      # clear any stale (discordant) RMC
        else:
            b["rmc"] = sorted(txs[best_tx], key=lambda x: x["s"])
            updated += 1
        payload = json.dumps(b, separators=(",", ":")).encode("utf-8")
        with gzip.GzipFile(path, "wb", mtime=0) as gz:
            gz.write(payload)
    print(f"  Injected RMC into {updated} bundles ({skipped_isoform} skipped: "
          f"isoform length mismatch)", flush=True)


if __name__ == "__main__":
    main()
