"""Build the UniProt -> {symbol, ensembl, chrom} cache from the hg38 file.

The all-19 aa_substitutions file is keyed only by UniProt (no gene symbol, no
genomic coordinates). We derive the symbol by taking each UniProt's genomic
locus from AlphaMissense_hg38.tsv.gz and overlapping it with MANE gene spans.
Result is cached to scripts/.uniprot_map.json so build_bundles can name the
UniProt-keyed data by gene symbol.

Usage:  python3 build_uniprot_map.py
"""
import json
import time

import config as C
import lib


def main():
    print("Loading MANE intervals…", flush=True)
    intervals = lib.load_mane_intervals(select_only=True)

    seen = {}  # uniprot -> (chrom, firstpos, ensembl)
    t0 = time.time()
    print(f"Scanning {C.ALPHAMISSENSE_TSV} for UniProt loci…", flush=True)
    with lib.openf(C.ALPHAMISSENSE_TSV) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 7)
            uniprot, transcript = f[5], f[6]
            if uniprot not in seen:
                seen[uniprot] = (f[0], int(f[1]), transcript)

    out = {}
    unmapped = 0
    for uniprot, (chrom, pos, ensembl) in seen.items():
        sym = lib.lookup_symbol(intervals, chrom, pos)
        if sym is None:
            unmapped += 1
            continue
        out[uniprot] = {"symbol": sym, "ensembl": ensembl, "chrom": chrom}

    with open(C.UNIPROT_MAP_JSON, "w") as fh:
        json.dump(out, fh, separators=(",", ":"))
    print(f"Wrote {C.UNIPROT_MAP_JSON}: {len(out)} UniProt→symbol "
          f"({unmapped} unmapped) in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
