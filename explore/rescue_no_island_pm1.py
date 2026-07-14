"""EXPLORATORY (not wired into the engine): can a secondary rule rescue no-island PM1
false negatives without introducing false positives?

Motivation: island-based PM1 is high-precision but misses ~575 no-island FN (dominated by
sub-threshold ClinVar hotspots, §3.3 of the paper). We ask whether any signal we already
compute — sub-threshold AlphaMissense, gnomAD regional missense constraint (RMC), cancer
hotspots, UniProt domains, or (as a circular baseline) ClinVar density — can, applied ONLY
to residues that overlap no island, flip enough true FN to TP without flooding TN into FP.

Data: engine PM1 calls + reasons from EREPO_SCORED_JSON; expert truth (PM1 in eRepo "Met")
from EREPO_TSV, matched by ClinVar Variation Id; per-residue signals from the shipped gene
bundles. Base benchmark: TP=618, FP=266, FN=604 (P 0.699, R 0.506, F1 0.587).

Headline result: no clean rescue exists. The no-island regime is ~7:1 TN:FN, and the FN we
most want (dense-ClinVar sub-threshold hotspots) live where the TN also cluster, so the
only signal that reaches them (ClinVar density) floods FP — which is itself evidence for
the island contiguity+length gate and the anti-circularity thesis. See rescue_findings.md.

Run:  python3 explore/rescue_no_island_pm1.py   (needs the same local data as scripts/)
"""
import csv
import gzip
import json
import os
import re
import statistics
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts"))
import config as C  # noqa: E402

PVAR = re.compile(r"\(p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\)")
AA3 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E",
    "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F",
    "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
}
TP, FP, FN_ALL = 618, 266, 604


def load_predictions():
    pred = {}
    with gzip.open(C.EREPO_SCORED_JSON, "rt") as fh:
        for line in fh:
            if not line.strip():
                continue
            d = json.loads(line)
            tags = d.get("link", {}).get("acmg_tags_v2", [])
            pm1 = next((t for t in tags if t.get("code") in ("PM1", "vPM1")), None)
            if pm1 is None:
                continue
            on = pm1.get("state") == "on"
            cids = pm1.get("comment_ids", []) or []
            reason = ("island" if on
                      else "notmiss" if any("NOT_MISSENSE" in c for c in cids)
                      else "noisland" if any("PM1_AM_OFF" in c for c in cids)
                      else "other")
            cvs = d.get("element", {}).get("extdb", {}).get("clinvar", []) or []
            if isinstance(cvs, (str, int)):
                cvs = [cvs]
            for c in cvs:
                pred[str(c)] = reason
    return pred


def load_candidates(pred):
    """No-island missense residues with expert truth + bundle-derived signals."""
    rows = []
    with open(C.EREPO_TSV) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            cid = (r.get("ClinVar Variation Id") or "").strip()
            if cid not in pred or pred[cid] != "noisland":
                continue
            m = PVAR.search(r.get("Variation", ""))
            if not m or m.group(1) not in AA3 or m.group(3) not in AA3:
                continue
            met = r.get("Applied Evidence Codes (Met)", "") or ""
            rows.append({"g": (r.get("HGNC Gene Symbol") or "").strip(),
                         "pos": int(m.group(2)), "ref": AA3[m.group(1)],
                         "truth": any(t.strip().startswith("PM1") for t in met.split(","))})
    return rows


def signals(rows):
    cache = {}

    def bundle(g):
        if g not in cache:
            p = os.path.join(C.GENES_DIR, g + ".json.gz")
            cache[g] = json.load(gzip.open(p, "rt")) if os.path.exists(p) else None
        return cache[g]

    out = []
    for r in rows:
        b = bundle(r["g"])
        if not b:
            continue
        pos, L = r["pos"], b["length"]
        if not (1 <= pos <= L) or b["ref"][pos - 1] != r["ref"]:
            continue                                   # isoform-concordance gate
        mean = b.get("mean") or []
        win = [mean[i] for i in range(max(0, pos - 6), min(L, pos + 5))
               if i < len(mean) and mean[i] is not None]
        isl = b.get("islands") or []
        oe = next((x["oe"] for x in (b.get("rmc") or [])
                   if x["s"] <= pos <= x["e"] and not x.get("whole")), None)
        r["s"] = {
            "am": mean[pos - 1] if pos - 1 < len(mean) and mean[pos - 1] is not None else 0,
            "am_wmean": statistics.mean(win) if win else 0,
            "idist": min((min(abs(pos - i["s"]), abs(pos - i["e"])) for i in isl), default=999),
            "oe": oe,
            "hot": any(abs(h["p"] - pos) <= 2 for h in (b.get("hotspots") or [])),
            "dom": any(d["s"] <= pos <= d["e"] for d in (b.get("domains") or [])),
            "pnear": sum(1 for v in b["clinvar"] if abs(v["p"] - pos) <= 10),
            "bnear": sum(1 for v in b["clinvar_benign"] if abs(v["p"] - pos) <= 10),
        }
        out.append(r)
    return out


def f1(p, r):
    return 2 * p * r / (p + r) if p + r else 0


def main():
    pred = load_predictions()
    cand = load_candidates(pred)
    mapped = signals(cand)
    fn = [r for r in mapped if r["truth"]]
    tn = [r for r in mapped if not r["truth"]]
    print(f"no-island candidates: {len(cand)} | mapped: {len(mapped)} "
          f"(FN={len(fn)}, TN={len(tn)}; ~{len(tn) / max(1, len(fn)):.0f}:1 TN:FN)")

    print("\n-- signal medians / rates (FN vs TN) --")
    for k in ("am", "am_wmean"):
        print(f"  {k:9}: FN {statistics.median(r['s'][k] for r in fn):.3f} | "
              f"TN {statistics.median(r['s'][k] for r in tn):.3f}")
    for k, lbl in (("oe", "RMC o/e (sub-regional)"), ):
        fv = [r['s'][k] for r in fn if r['s'][k] is not None]
        tv = [r['s'][k] for r in tn if r['s'][k] is not None]
        print(f"  {lbl}: FN {statistics.median(fv):.3f} | TN {statistics.median(tv):.3f}")
    for k in ("hot", "dom"):
        print(f"  {k:9}: FN {sum(r['s'][k] for r in fn)}/{len(fn)} | "
              f"TN {sum(r['s'][k] for r in tn)}/{len(tn)}")

    def rule(name, f):
        rf = sum(1 for r in fn if f(r["s"]))
        nf = sum(1 for r in tn if f(r["s"]))
        p, rec = (TP + rf) / (TP + rf + FP + nf), (TP + rf) / (TP + FN_ALL)
        print(f"  {name:48} rescue={rf:3d} newFP={nf:4d}  P={p:.3f} R={rec:.3f} F1={f1(p, rec):.3f}")

    print(f"\n-- rescue rules (base P=0.699 R=0.506 F1=0.587) --")
    rule("raw AM residue >0.5", lambda s: s['am'] > 0.5)
    rule("UniProt domain member", lambda s: s['dom'])
    rule("[circular] >=3 P/LP & 0 benign +-10", lambda s: s['pnear'] >= 3 and s['bnear'] == 0)
    rule("gnomAD RMC o/e <0.6", lambda s: s['oe'] is not None and s['oe'] < 0.6)
    rule("just-outside-island <=8aa", lambda s: s['idist'] <= 8)
    rule("cancer-hotspot +-2aa", lambda s: s['hot'])
    rule("AM wmean>0.55 AND (hotspot or RMC<0.5)",
         lambda s: s['am_wmean'] > 0.55 and (s['hot'] or (s['oe'] is not None and s['oe'] < 0.5)))
    rule("stacked safe (hotspot | outside<=8 | AM&RMC)",
         lambda s: s['hot'] or s['idist'] <= 8
         or (s['am_wmean'] > 0.55 and s['oe'] is not None and s['oe'] < 0.5))

    sub = [r for r in fn if r["s"]["pnear"] >= 3]
    print(f"\nsub-threshold-hotspot-like FN (>=3 P/LP nearby): {len(sub)}/{len(fn)}; reachable by "
          f"cancer-hotspot={sum(r['s']['hot'] for r in sub)}, "
          f"near-island<=8={sum(r['s']['idist'] <= 8 for r in sub)}, "
          f"RMC<0.6={sum(r['s']['oe'] is not None and r['s']['oe'] < 0.6 for r in sub)}")


if __name__ == "__main__":
    main()
