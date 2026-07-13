"""Stage 3 — emit public/data/index.json from the gene bundles (drives search)."""
import glob
import gzip
import json
import os

import config as C


def main():
    files = sorted(glob.glob(os.path.join(C.GENES_DIR, "*.json.gz")))
    index = []
    for path in files:
        with gzip.open(path, "rt") as fh:
            b = json.load(fh)
        index.append({"symbol": b["gene"], "uniprot": b["uniprot"], "len": b["length"]})
    index.sort(key=lambda d: d["symbol"])
    os.makedirs(C.DATA_DIR, exist_ok=True)
    with open(C.INDEX_JSON, "w") as fh:
        json.dump(index, fh, separators=(",", ":"))
    print(f"Wrote {C.INDEX_JSON} with {len(index)} genes", flush=True)


if __name__ == "__main__":
    main()
