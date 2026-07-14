"""Build the PM1 benchmark dataset by replicating the SeqOne eRepo benchmark.

Mirrors benchmark_ACMG/benchmark_acmg.py exactly:
  - prediction = the engine's acmg_tags_v2 PM1 state (island-based, scored with
    --am-islands) from erepo_variants_updated.json.gz,
  - ground truth = PM1 (strength-independent) in the eRepo "Applied Evidence
    Codes (Met)",
  - matched by ClinVar Variation Id,
  - TP/FP/FN over the matched set.

This reproduces the headline PM1 numbers (the engine uses genomic island overlap,
so there is no isoform mapping and no benign downgrade). Output:
public/data/pm1_benchmark.json.

Usage:  python3 build_pm1_benchmark.py
"""
import csv
import gzip
import json
import os
import re
import sys

import config as C

PVAR_RE = re.compile(r"\(p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter|=)\)")
AA3 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E",
    "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F",
    "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
}


def parse_predictions(path):
    """clinvar_id -> {'pred': bool, 'pstr': str} from the engine's PM1 tag."""
    pred = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if not line.strip():
                continue
            d = json.loads(line)
            tags = d.get("link", {}).get("acmg_tags_v2", [])
            pm1 = next((t for t in tags if t.get("code") in ("PM1", "vPM1")), None)
            if pm1 is None:
                continue
            on = pm1.get("state") == "on"
            pstr = (pm1.get("evidence_level") or "").title() if on else ""
            cvids = d.get("element", {}).get("extdb", {}).get("clinvar", []) or []
            if isinstance(cvids, (str, int)):
                cvids = [cvids]
            for cid in cvids:
                pred[str(cid)] = {"pred": on, "pstr": pstr}
    return pred


def clean_variation(var):
    """A short human label from the eRepo Variation string."""
    m = re.search(r"\(([^)]+)\):(c\.[^ ]+)(?: \((p\.[^)]+)\))?", var)
    if m:
        return f"{m.group(2)}" + (f" ({m.group(3)})" if m.group(3) else "")
    return var[:60]


def main():
    for p in (C.EREPO_TSV, C.EREPO_SCORED_JSON):
        if not os.path.exists(p):
            print(f"missing: {p}", file=sys.stderr); sys.exit(1)

    print(f"Parsing engine predictions {C.EREPO_SCORED_JSON}…", flush=True)
    pred = parse_predictions(C.EREPO_SCORED_JSON)
    print(f"  {len(pred)} ClinVar-id predictions", flush=True)

    rows, counts = [], {"tp": 0, "fp": 0, "fn": 0, "tn": 0}
    unmatched = 0
    with open(C.EREPO_TSV) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            cid = (r.get("ClinVar Variation Id") or "").strip()
            if not cid or cid not in pred:
                unmatched += 1
                continue
            met = r.get("Applied Evidence Codes (Met)", "")
            met_list = [t.strip() for t in met.split(",") if t.strip()]
            truth = any(t.startswith("PM1") for t in met_list)
            tstr = next((t for t in met_list if t.startswith("PM1")), "")
            p = pred[cid]
            cat = ("tp" if truth and p["pred"] else "fp" if p["pred"] and not truth
                   else "fn" if truth and not p["pred"] else "tn")
            counts[cat] += 1
            if cat == "tn":
                continue

            gene = (r.get("HGNC Gene Symbol") or "").strip()
            variation = r.get("Variation", "")
            m = PVAR_RE.search(variation)
            vshort = f"{AA3.get(m.group(1), '?')}{m.group(2)}{AA3.get(m.group(3), m.group(3))}" if m else ""
            rows.append({
                "g": gene,
                "v": vshort,                       # 1-letter p. (for the app PM1 link), may be ""
                "hgvs": clean_variation(variation),
                "dis": (r.get("Disease") or "").strip()[:90],
                "panel": (r.get("Expert Panel") or "").strip(),
                "link": (r.get("Evidence Repo Link") or "").strip(),
                "clinvar": cid,
                "assertion": (r.get("Assertion") or "").strip(),
                "met": ", ".join(met_list),
                "truth": truth, "tstr": tstr,
                "pred": p["pred"], "pstr": p["pstr"],
                "cat": cat,
            })

    tp, fp, fn = counts["tp"], counts["fp"], counts["fn"]
    prec = tp / (tp + fp) if tp + fp else 0
    rec = tp / (tp + fn) if tp + fn else 0
    f1 = 2 * prec * rec / (prec + rec) if prec + rec else 0
    summary = {**counts, "precision": round(prec, 4), "recall": round(rec, 4),
               "f1": round(f1, 4), "matched": tp + fp + fn + counts["tn"],
               "n_rows": len(rows)}
    rows.sort(key=lambda x: (x["cat"], x["g"], x["hgvs"]))
    os.makedirs(C.DATA_DIR, exist_ok=True)
    with open(C.PM1_BENCHMARK_JSON, "w") as fh:
        json.dump({"summary": summary, "variants": rows}, fh, separators=(",", ":"))
    print("PM1 benchmark:", summary, flush=True)
    print(f"Wrote {C.PM1_BENCHMARK_JSON}", flush=True)


if __name__ == "__main__":
    main()
