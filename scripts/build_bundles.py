"""Stage 1 — build per-gene bundles from AlphaMissense all-19 substitutions.

Streams AlphaMissense_aa_substitutions.tsv.gz (keyed by UniProt, all 19
substitutions per residue → dense heatmap). Groups by UniProt (the file is
sorted by uniprot_id) and flushes each protein as its own bundle. Gene symbol,
Ensembl transcript come from the UniProt map cache (build_uniprot_map.py);
islands + domains are joined by UniProt. ClinVar is added by add_clinvar.py.

Usage:
  python3 build_bundles.py [--genes KIF1A,TP53] [--aasub PATH]
"""
import argparse
import gzip
import json
import os
import sys
import time
from collections import defaultdict

import config as C
import lib


def write_bundle(uniprot, acc, umap, islands_u, domains, targets, stats):
    rec = umap.get(uniprot)
    if rec is None:
        stats["nosymbol"] += 1
        return
    symbol = rec["symbol"]
    if targets is not None and symbol not in targets:
        return
    ref_by_pos, scores_by_pos = acc["ref"], acc["scores"]
    length = acc["maxpos"]
    ref_str = "".join(ref_by_pos.get(p, "X") for p in range(1, length + 1))
    mean = [
        (round(sum(s.values()) / len(s), 4) if (s := scores_by_pos.get(p)) else None)
        for p in range(1, length + 1)
    ]
    bundle = {
        "gene": symbol,
        "uniprot": uniprot,
        "ensembl": rec.get("ensembl", ""),
        "refseq": "",
        "length": length,
        "aa_order": C.AA_ORDER,
        "mean": mean,
        "ref": ref_str,
        "grid": lib.pack_grid(length, ref_by_pos, scores_by_pos),
        "islands": islands_u.get(uniprot, []),
        "domains": domains.get(uniprot, []),
        "clinvar": [],
        "clinvar_benign": [],
    }
    path = os.path.join(C.GENES_DIR, symbol + ".json.gz")
    payload = json.dumps(bundle, separators=(",", ":")).encode("utf-8")
    with gzip.GzipFile(path, "wb", mtime=0) as gz:
        gz.write(payload)
    stats["written"] += 1


def find_columns(line):
    fields = [x.strip() for x in line.lstrip("#").rstrip("\n").split("\t")]
    col = {n: i for i, n in enumerate(fields)}
    return col["uniprot_id"], col["protein_variant"], col["am_pathogenicity"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes")
    ap.add_argument("--aasub", default=C.ALPHAMISSENSE_AASUB)
    args = ap.parse_args()
    targets = set(args.genes.split(",")) if args.genes else None

    if not os.path.exists(C.UNIPROT_MAP_JSON):
        print("Missing UniProt map — run build_uniprot_map.py first.", file=sys.stderr)
        sys.exit(1)

    os.makedirs(C.GENES_DIR, exist_ok=True)
    print("Loading UniProt map / islands / domains…", flush=True)
    umap = json.load(open(C.UNIPROT_MAP_JSON))
    islands_u = lib.load_islands_by_uniprot()
    domains = lib.load_domains()
    print(f"  UniProt map: {len(umap)} | island proteins: {len(islands_u)} | "
          f"domain proteins: {len(domains)}", flush=True)

    stats = {"written": 0, "nosymbol": 0, "rows": 0}
    t0 = time.time()
    i_u = i_pv = i_sc = None
    cur = None
    acc = None
    flushed = set()

    print(f"Streaming {args.aasub}…", flush=True)
    with lib.openf(args.aasub) as fh:
        for line in fh:
            if line.startswith("#"):
                if i_u is None and "protein_variant" in line:
                    i_u, i_pv, i_sc = find_columns(line)
                continue
            if i_u is None:
                if "protein_variant" in line:  # header without leading '#'
                    i_u, i_pv, i_sc = find_columns(line)
                    continue
                raise RuntimeError("aa_substitutions header not found")

            f = line.rstrip("\n").split("\t")
            uniprot = f[i_u]
            if uniprot != cur:
                if cur is not None:
                    write_bundle(cur, acc, umap, islands_u, domains, targets, stats)
                    flushed.add(cur)
                    if stats["written"] % 2000 == 0 and stats["written"]:
                        print(f"  {stats['written']} bundles ({time.time()-t0:.0f}s)", flush=True)
                if uniprot in flushed:
                    print(f"  WARNING: {uniprot} reappeared (file not uniprot-sorted)", file=sys.stderr)
                cur = uniprot
                acc = {"ref": {}, "scores": defaultdict(dict), "maxpos": 0}

            parsed = lib.parse_pvar(f[i_pv])
            if parsed is None:
                continue
            ref, pos, alt = parsed
            if ref not in C.AA_INDEX or alt not in C.AA_INDEX:
                continue
            stats["rows"] += 1
            acc["ref"][pos] = ref
            acc["scores"][pos][alt] = float(f[i_sc])
            if pos > acc["maxpos"]:
                acc["maxpos"] = pos

    if cur is not None:
        write_bundle(cur, acc, umap, islands_u, domains, targets, stats)

    print(f"Done: {stats['written']} bundles, {stats['rows']} substitutions, "
          f"{stats['nosymbol']} UniProt without a mapped symbol, "
          f"{time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
