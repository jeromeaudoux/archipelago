"""Derive per-transcript exon boundaries (in residue space) and inject into bundles.

Uses AlphaMissense_hg38.tsv.gz (genomic position per residue, same isoform as the
bundles): within a coding exon consecutive residues' codons step by 3 nt; a larger
jump marks an intron, i.e. the start of a new exon. Stores `exons` = list of
exon-start residues (1-based) for multi-exon transcripts.

Usage:  python3 build_exons.py
"""
import gzip
import json
import os
import sys
import time
from collections import defaultdict

import config as C
import lib

STEP_MAX = 5  # within-exon codon step is 3 nt; anything larger is an intron boundary


def flush(genes, umap, stats):
    for uniprot, minpos in genes.items():
        rec = umap.get(uniprot)
        if rec is None:
            continue
        path = os.path.join(C.GENES_DIR, rec["symbol"] + ".json.gz")
        if not os.path.exists(path):
            continue
        res = sorted(minpos)
        starts = [res[0]]
        for i in range(1, len(res)):
            if abs(minpos[res[i]] - minpos[res[i - 1]]) > STEP_MAX or res[i] - res[i - 1] > 1:
                starts.append(res[i])
        with gzip.open(path, "rt") as fh:
            b = json.load(fh)
        if len(starts) > 1 and starts[-1] <= b["length"]:
            b["exons"] = starts
        elif "exons" in b:
            del b["exons"]
        else:
            continue
        payload = json.dumps(b, separators=(",", ":")).encode("utf-8")
        with gzip.GzipFile(path, "wb", mtime=0) as gz:
            gz.write(payload)
        stats["written"] += 1


def main():
    if not os.path.exists(C.UNIPROT_MAP_JSON):
        print("Missing UniProt map — run build_uniprot_map.py first.", file=sys.stderr); sys.exit(1)
    umap = json.load(open(C.UNIPROT_MAP_JSON))
    stats = {"written": 0}
    t0 = time.time()
    cur_chrom = None
    genes = defaultdict(lambda: defaultdict(lambda: 1 << 62))  # uniprot -> residue -> min POS

    print(f"Streaming {C.ALPHAMISSENSE_TSV} for exon boundaries…", flush=True)
    with lib.openf(C.ALPHAMISSENSE_TSV) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 8)
            chrom, pos, uniprot, transcript, pvar = f[0], int(f[1]), f[5], f[6], f[7]
            if uniprot not in umap:
                continue
            if chrom != cur_chrom:
                if cur_chrom is not None:
                    flush(genes, umap, stats)
                    genes = defaultdict(lambda: defaultdict(lambda: 1 << 62))
                    print(f"  {cur_chrom}: {stats['written']} genes ({time.time()-t0:.0f}s)", flush=True)
                cur_chrom = chrom
            parsed = lib.parse_pvar(pvar)
            if parsed is None:
                continue
            r = parsed[1]
            d = genes[uniprot]
            if pos < d[r]:
                d[r] = pos
    if cur_chrom is not None:
        flush(genes, umap, stats)
    print(f"Done: exons injected into {stats['written']} bundles ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
