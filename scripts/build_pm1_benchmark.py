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
from collections import Counter

import config as C

PVAR_RE = re.compile(r"\(p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter|=)\)")
VLETTER_RE = re.compile(r"^([A-Z])(\d+)([A-Z*=])$")
CNUM_RE = re.compile(r"c\.(\d+)")

# Window (± residues) used to characterise the ClinVar neighbourhood of a variant.
WIN = 10


def load_bundle(gene, cache):
    if gene not in cache:
        path = os.path.join(C.GENES_DIR, gene + ".json.gz")
        b = None
        if os.path.exists(path):
            with gzip.open(path, "rt") as fh:
                b = json.load(fh)
        cache[gene] = b
    return cache[gene]


def variant_pos(row):
    """(residue, ref_letter|None) from the 1-letter p. or, failing that, the c. HGVS."""
    m = VLETTER_RE.match(row["v"])
    if m:
        return int(m.group(2)), m.group(1)
    cm = CNUM_RE.search(row["c"])
    if cm:
        return (int(cm.group(1)) + 2) // 3, None
    return None, None


def enrich(row, cache):
    """Attach an honest divergence category + supporting signals to an FP/FN row.

    Signals come from the gene's AlphaMissense bundle (AM islands, per-residue
    ClinVar P/LP and B/LB, mean AM, domains). The eRepo variant may sit on a
    different isoform than the AM canonical bundle; when the reference amino acid
    disagrees we mark the row `unmapped` rather than guess.
    """
    if row["cat"] == "fn" and row["pr"] == "notmiss":
        row["dcat"] = "fn_not_missense"       # engine-authoritative, no bundle needed
        return
    b = load_bundle(row["g"], cache)
    pos, refaa = variant_pos(row)
    if b is None or pos is None or not (1 <= pos <= b["length"]):
        row["dcat"] = "unmapped"
        return
    if refaa is not None and b["ref"][pos - 1] != refaa:
        row["dcat"] = "unmapped"
        return

    islands = b.get("islands", [])
    in_isl = any(i["s"] <= pos <= i["e"] for i in islands)
    idist = 0 if in_isl else min(
        (min(abs(pos - i["s"]), abs(pos - i["e"])) for i in islands), default=-1)
    pnear = sum(1 for r in b["clinvar"] if abs(r["p"] - pos) <= WIN)
    bnear = sum(1 for r in b["clinvar_benign"] if abs(r["p"] - pos) <= WIN)
    mean = b.get("mean") or []
    am = mean[pos - 1] if pos - 1 < len(mean) and mean[pos - 1] is not None else None
    dom = any(d["s"] <= pos <= d["e"] for d in b.get("domains", []))
    row.update({"pnear": pnear, "bnear": bnear, "idist": idist, "dom": dom,
                "amv": round(am, 2) if am is not None else None})

    if row["cat"] == "fp":
        if bnear >= 2:
            row["dcat"] = "fp_benign_conflict"
        elif pnear >= 3:
            row["dcat"] = "fp_clinvar_hotspot"
        else:
            row["dcat"] = "fp_am_only"
    else:  # fn
        if pnear >= 3:
            row["dcat"] = "fn_subthreshold_hotspot"
        elif 1 <= idist <= WIN:
            row["dcat"] = "fn_near_island"
        else:
            row["dcat"] = "fn_weak_am"
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
            cids = pm1.get("comment_ids", []) or []
            if on:
                reason = "island"
            elif any("NOT_MISSENSE" in c for c in cids):
                reason = "notmiss"     # excluded: not a missense / in-frame variant
            elif any("PM1_AM_OFF" in c for c in cids):
                reason = "noisland"    # missense but no island overlap
            else:
                reason = "other"
            cvids = d.get("element", {}).get("extdb", {}).get("clinvar", []) or []
            if isinstance(cvids, (str, int)):
                cvids = [cvids]
            for cid in cvids:
                pred[str(cid)] = {"pred": on, "pstr": pstr, "pr": reason}
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
            cm = re.search(r":(c\.[^ )]+)", variation) or re.search(r"(c\.[0-9][^ ),]*)", r.get("HGVS Expressions", ""))
            rows.append({
                "g": gene,
                "v": vshort,                       # 1-letter p. (for the app PM1 link), may be ""
                "c": cm.group(1) if cm else "",    # c. HGVS (fallback for the PM1 view)
                "hgvs": clean_variation(variation),
                "dis": (r.get("Disease") or "").strip()[:90],
                "panel": (r.get("Expert Panel") or "").strip(),
                "link": (r.get("Evidence Repo Link") or "").strip(),
                "clinvar": cid,
                "assertion": (r.get("Assertion") or "").strip(),
                "met": ", ".join(met_list),
                "truth": truth, "tstr": tstr,
                "pred": p["pred"], "pstr": p["pstr"], "pr": p.get("pr", "other"),
                "cat": cat,
            })

    bundle_cache = {}
    for r in rows:
        if r["cat"] in ("fp", "fn"):
            enrich(r, bundle_cache)

    tp, fp, fn = counts["tp"], counts["fp"], counts["fn"]
    prec = tp / (tp + fp) if tp + fp else 0
    rec = tp / (tp + fn) if tp + fn else 0
    f1 = 2 * prec * rec / (prec + rec) if prec + rec else 0
    fn_notmiss = sum(1 for r in rows if r["cat"] == "fn" and r["pr"] == "notmiss")
    fn_noisland = sum(1 for r in rows if r["cat"] == "fn" and r["pr"] == "noisland")
    dcats = Counter(r["dcat"] for r in rows if r["cat"] in ("fp", "fn"))
    summary = {**counts, "precision": round(prec, 4), "recall": round(rec, 4),
               "f1": round(f1, 4), "matched": tp + fp + fn + counts["tn"],
               "n_rows": len(rows),
               "fn_not_missense": fn_notmiss, "fn_no_island": fn_noisland,
               "fn_other": fn - fn_notmiss - fn_noisland,
               "dcats": dict(dcats)}
    rows.sort(key=lambda x: (x["cat"], x["g"], x["hgvs"]))
    os.makedirs(C.DATA_DIR, exist_ok=True)
    with open(C.PM1_BENCHMARK_JSON, "w") as fh:
        json.dump({"summary": summary, "variants": rows}, fh, separators=(",", ":"))
    print("PM1 benchmark:", summary, flush=True)
    print(f"Wrote {C.PM1_BENCHMARK_JSON}", flush=True)


if __name__ == "__main__":
    main()
