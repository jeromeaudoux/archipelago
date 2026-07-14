"""Inject gnomAD regional missense constraint (RMC) into gene bundles.

Two gnomAD v2.1.1 inputs:
  - transcripts_with_rmc:    per-transcript sub-regions with observed/expected.
  - transcripts_without_rmc: one gene-level o/e (no significant sub-region).

Gene identity is resolved via Ensembl gene id (gene_id) → MANE Select symbol
(robust to symbol drift), falling back to the file's gene_name. Sub-regions are
placed in residue space on the transcript whose length matches the AlphaMissense
isoform (or whose Ensembl id matches it); a gene with only a discordant isoform is
skipped. Gene-level o/e (without-RMC) is drawn as a single whole-protein band.

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
import lib

RESI_RE = re.compile(r"(\d+)$")


def load_mane_ensg():
    """Ensembl gene id (base) -> MANE Select symbol."""
    m = {}
    with lib.openf(C.MANE_SUMMARY) as fh:
        hdr = fh.readline().lstrip("#").rstrip("\n").split("\t")
        col = {n: i for i, n in enumerate(hdr)}
        for line in fh:
            v = line.rstrip("\n").split("\t")
            if v[col["MANE_status"]] != "MANE Select":
                continue
            m[lib.base(v[col["Ensembl_Gene"]])] = v[col["symbol"]]
    return m


def load_with_rmc(sym_of):
    """symbol -> transcript_base -> list of regions."""
    by_gene = defaultdict(lambda: defaultdict(list))
    with open(C.GNOMAD_RMC_TSV) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            sym = sym_of(r["gene_id"], r["gene_name"])
            ms, me = RESI_RE.search(r["start_aa"] or ""), RESI_RE.search(r["stop_aa"] or "")
            if not (sym and ms and me):
                continue
            try:
                oe = float(r["oe"]); obs = int(float(r["obs"])); exp = float(r["exp"]); pv = float(r["p_value"])
            except (ValueError, KeyError):
                continue
            a, c = int(ms.group(1)), int(me.group(1))
            by_gene[sym][lib.base(r["transcript"])].append(
                {"s": min(a, c), "e": max(a, c), "oe": round(oe, 4), "obs": obs, "exp": round(exp, 1), "p": pv})
    return by_gene


def load_without_rmc(sym_of):
    """symbol -> {oe, obs, exp} gene-level constraint."""
    out = {}
    with open(C.GNOMAD_NO_RMC_TSV) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            sym = sym_of(r["gene_id"], r["gene_name"])
            if not sym or sym in out:
                continue
            try:
                out[sym] = {"oe": round(float(r["oe"]), 4), "obs": int(float(r["obs"])), "exp": round(float(r["exp"]), 1)}
            except (ValueError, KeyError):
                continue
    return out


def pick_transcript(txs, ensembl_base, length):
    """Choose the RMC transcript matching the AM isoform (id match, else closest length)."""
    if ensembl_base in txs:
        return ensembl_base
    best, best_fit = None, None
    for tx, regs in txs.items():
        maxaa = max(r["e"] for r in regs)
        if maxaa <= length + 1 and (best_fit is None or abs(maxaa - length) < best_fit):
            best, best_fit = tx, abs(maxaa - length)
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes")
    args = ap.parse_args()
    targets = set(args.genes.split(",")) if args.genes else None

    if not os.path.isdir(C.GENES_DIR):
        print("No bundles found — run build_bundles.py first.", file=sys.stderr); sys.exit(1)

    mane_ensg = load_mane_ensg()

    def sym_of(gene_id, gene_name):
        return mane_ensg.get(lib.base(gene_id)) or (gene_name or "").strip()

    print("Reading gnomAD RMC (with + without)…", flush=True)
    with_rmc = load_with_rmc(sym_of)
    gene_oe = load_without_rmc(sym_of)
    print(f"  with-RMC genes: {len(with_rmc)} | gene-level genes: {len(gene_oe)}", flush=True)

    syms = set(with_rmc) | set(gene_oe)
    updated = level = skipped = cleared = 0
    for sym in syms:
        if targets is not None and sym not in targets:
            continue
        path = os.path.join(C.GENES_DIR, sym + ".json.gz")
        if not os.path.exists(path):
            continue
        with gzip.open(path, "rt") as fh:
            b = json.load(fh)
        length = b["length"]; enst = lib.base(b.get("ensembl", ""))
        rmc = None
        if sym in with_rmc:
            tx = pick_transcript(with_rmc[sym], enst, length)
            if tx is not None:
                rmc = sorted(with_rmc[sym][tx], key=lambda x: x["s"])
                updated += 1
        if rmc is None and sym in gene_oe:            # gene-level fallback
            g = gene_oe[sym]
            rmc = [{"s": 1, "e": length, "oe": g["oe"], "obs": g["obs"], "exp": g["exp"], "whole": 1}]
            level += 1
        if rmc is None:
            skipped += 1
            if "rmc" not in b:
                continue
            del b["rmc"]; cleared += 1
        else:
            b["rmc"] = rmc
        payload = json.dumps(b, separators=(",", ":")).encode("utf-8")
        with gzip.GzipFile(path, "wb", mtime=0) as gz:
            gz.write(payload)
    print(f"  sub-regional: {updated} | gene-level: {level} | skipped (isoform): {skipped} "
          f"| stale cleared: {cleared}", flush=True)


if __name__ == "__main__":
    main()
