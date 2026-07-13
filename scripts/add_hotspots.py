"""Inject cancerhotspots.org residue hotspots into gene bundles.

Reads the "Hotspot_Residues" sheet of hotspots_v3.xlsx and, for genes that have
entries, adds a `hotspots` array to the bundle (dropped for genes without any).
Positions are validated against the AlphaMissense reference sequence.

Usage:  python3 add_hotspots.py [--genes TP53,BRAF]
"""
import argparse
import glob
import gzip
import json
import os
import re
import sys
from collections import defaultdict

import config as C

CODON_RE = re.compile(r"^([A-Za-z])(\d+)$")


def load_hotspots():
    import openpyxl
    wb = openpyxl.load_workbook(C.CANCERHOTSPOTS_XLSX, read_only=True, data_only=True)
    ws = wb["Hotspot_Residues"]
    rows = ws.iter_rows(values_only=True)
    hdr = list(next(rows))
    col = {name: i for i, name in enumerate(hdr) if name}
    ci_g = col["Hugo_Symbol"]; ci_c = col["Codon"]; ci_p = col["Codon_Position"]
    ci_q = col.get("Q value"); ci_onc = col.get("In OncoKB v3.14"); ci_n = col.get("# mut in MSK")

    by_gene = defaultdict(list)
    for r in rows:
        sym = r[ci_g]
        if not sym:
            continue
        pos = r[ci_p]
        ref = ""
        m = CODON_RE.match(str(r[ci_c]).strip()) if r[ci_c] else None
        if m:
            ref = m.group(1)
            if pos is None:
                pos = int(m.group(2))
        if pos is None:
            continue
        rec = {"p": int(pos), "r": ref}
        if ci_q is not None and r[ci_q] is not None:
            rec["q"] = float(r[ci_q])
        if ci_n is not None and r[ci_n] is not None:
            rec["n"] = int(r[ci_n])
        rec["onc"] = 1 if (ci_onc is not None and r[ci_onc]) else 0
        by_gene[sym].append(rec)
    for v in by_gene.values():
        v.sort(key=lambda d: d["p"])
    return by_gene


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes")
    args = ap.parse_args()
    targets = set(args.genes.split(",")) if args.genes else None

    if not os.path.isdir(C.GENES_DIR):
        print("No bundles found — run build_bundles.py first.", file=sys.stderr)
        sys.exit(1)
    print(f"Reading {C.CANCERHOTSPOTS_XLSX}…", flush=True)
    by_gene = load_hotspots()
    print(f"  {sum(len(v) for v in by_gene.values())} hotspot residues across "
          f"{len(by_gene)} genes", flush=True)

    updated = mism = 0
    for sym, recs in by_gene.items():
        if targets is not None and sym not in targets:
            continue
        path = os.path.join(C.GENES_DIR, sym + ".json.gz")
        if not os.path.exists(path):
            continue
        with gzip.open(path, "rt") as fh:
            b = json.load(fh)
        ref_str, length = b["ref"], b["length"]
        kept = []
        for rec in recs:
            p = rec["p"]
            if p < 1 or p > length:
                mism += 1; continue
            if rec["r"] and ref_str[p - 1] != rec["r"]:
                mism += 1; continue          # off-isoform position
            kept.append(rec)
        if not kept:
            continue
        b["hotspots"] = kept
        payload = json.dumps(b, separators=(",", ":")).encode("utf-8")
        with gzip.GzipFile(path, "wb", mtime=0) as gz:
            gz.write(payload)
        updated += 1
    print(f"  Injected cancer hotspots into {updated} gene bundles "
          f"({mism} off-isoform positions dropped)", flush=True)


if __name__ == "__main__":
    main()
