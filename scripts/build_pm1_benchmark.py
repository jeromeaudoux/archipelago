"""Build the PM1 benchmark dataset: ClinGen-curated PM1 (eRepo) vs Archipelago
island-PM1, per variant, classified TP/FP/FN.

Ground truth: a variant is PM1-positive if any "Applied Evidence Codes (Met)" token
starts with PM1. Prediction: the residue falls in an AlphaMissense island (and the
island is not benign-contradicted), matching the app's PM1 rule. Only missense
variants that map onto an AlphaMissense isoform (ref aa concordant) are scored.

Output: public/data/pm1_benchmark.json (TP/FP/FN rows + summary counts).

Usage:  python3 build_pm1_benchmark.py
"""
import csv
import gzip
import json
import os
import re
import sys

import config as C

AA3 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E",
    "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F",
    "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V", "Ter": "*",
}
PVAR_RE = re.compile(r"\(p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter|=)\)")

_bundle_cache: dict = {}


def load_gene(symbol):
    if symbol in _bundle_cache:
        return _bundle_cache[symbol]
    path = os.path.join(C.GENES_DIR, symbol + ".json.gz")
    b = None
    if os.path.exists(path):
        with gzip.open(path, "rt") as fh:
            full = json.load(fh)
        b = {"ref": full["ref"], "length": full["length"], "islands": full["islands"]}
    _bundle_cache[symbol] = b
    return b


def predict(b, pos):
    """Return (pred_bool, strength_label, island[s,e] | None) mirroring pm1.ts."""
    isl = next((i for i in b["islands"] if i["s"] <= pos <= i["e"]), None)
    if not isl:
        return False, "Not met", None
    plp, blb = isl.get("plp", 0), isl.get("blb", 0)
    span = [isl["s"], isl["e"]]
    if blb > 0 and blb >= plp:
        return False, "Not met (benign-heavy)", span
    if plp >= 10 and blb == 0:
        return True, "Strong", span
    if plp >= 3 and blb == 0:
        return True, "Moderate", span
    return True, "Supporting", span


def main():
    if not os.path.exists(C.EREPO_TSV):
        print(f"eRepo not found: {C.EREPO_TSV}", file=sys.stderr); sys.exit(1)

    rows = []
    counts = {"tp": 0, "fp": 0, "fn": 0, "tn": 0}
    skipped_nogene = skipped_nopvar = refmismatch = 0

    with open(C.EREPO_TSV) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            gene = (r.get("HGNC Gene Symbol") or "").strip()
            met = r.get("Applied Evidence Codes (Met)", "")
            truth = any(t.strip().startswith("PM1") for t in met.split(","))
            truth_str = next((t.strip() for t in met.split(",") if t.strip().startswith("PM1")), "")

            m = PVAR_RE.search(r.get("Variation", ""))
            if not m:
                skipped_nopvar += 1; continue
            ref, alt = AA3.get(m.group(1)), (m.group(3) if m.group(3) in ("=", "Ter") else AA3.get(m.group(3)))
            pos = int(m.group(2))
            if not ref or alt in ("*", "=", None):   # missense only
                skipped_nopvar += 1; continue

            b = load_gene(gene)
            if b is None:
                skipped_nogene += 1; continue
            if pos < 1 or pos > b["length"] or b["ref"][pos - 1] != ref:
                refmismatch += 1; continue

            pred, pred_str, span = predict(b, pos)
            cat = ("tp" if truth and pred else "fp" if pred and not truth
                   else "fn" if truth and not pred else "tn")
            counts[cat] += 1
            if cat == "tn":
                continue  # too many & uninformative; summarized only

            rows.append({
                "g": gene, "v": f"{ref}{pos}{alt}", "p": pos,
                "dis": (r.get("Disease") or "").strip()[:80],
                "panel": (r.get("Expert Panel") or "").strip(),
                "link": (r.get("Evidence Repo Link") or "").strip(),
                "clinvar": (r.get("ClinVar Variation Id") or "").strip(),
                "truth": truth, "tstr": truth_str,
                "pred": pred, "pstr": pred_str,
                "isl": span, "cat": cat,
            })

    tp, fp, fn = counts["tp"], counts["fp"], counts["fn"]
    prec = tp / (tp + fp) if tp + fp else 0
    rec = tp / (tp + fn) if tp + fn else 0
    f1 = 2 * prec * rec / (prec + rec) if prec + rec else 0
    summary = {**counts, "precision": round(prec, 4), "recall": round(rec, 4),
               "f1": round(f1, 4), "n_rows": len(rows),
               "skipped_no_gene": skipped_nogene, "skipped_non_missense": skipped_nopvar,
               "ref_mismatch": refmismatch}
    rows.sort(key=lambda x: (x["cat"], x["g"], x["p"]))
    os.makedirs(C.DATA_DIR, exist_ok=True)
    with open(C.PM1_BENCHMARK_JSON, "w") as fh:
        json.dump({"summary": summary, "variants": rows}, fh, separators=(",", ":"))
    print("PM1 benchmark:", summary, flush=True)
    print(f"Wrote {C.PM1_BENCHMARK_JSON}", flush=True)


if __name__ == "__main__":
    main()
