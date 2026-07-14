"""Quantify how (in)consistently ClinGen expert panels apply PM1 to neighbouring variants.

Motivation: our island-based PM1 is benchmarked against ClinGen eRepo expert curation as
"ground truth", which implicitly assumes that PM1 label is consistent. It is not. Even a
single Variant Curation Expert Panel (VCEP), working from a refined gene-specific ACMG
specification, applies PM1 to one missense variant and withholds it from a near-identical
neighbour.

We scan eRepo missense variants and, within the SAME gene and SAME expert panel, look for
pairs where one variant has PM1 applied and a nearby variant does not — reported in two
tiers of rigour:

  Tier 1 (airtight): both variants Pathogenic/Likely-Pathogenic (so PM1 applicability is
    comparable) and the negative one has PM1 *explicitly* in "Applied Evidence Codes
    (Not Met)" (so "they didn't evaluate PM1" cannot explain it).
  Tier 2 (prevalence): the negative neighbour merely lacks PM1 (absent, looser).

Also flags any single variant curated more than once by the same panel with discordant
PM1 — the strongest possible inconsistency.

Outputs public/data/pm1_reproducibility.json and paper/supp_pm1_reproducibility.md.

Usage:  python3 analyze_pm1_reproducibility.py
"""
import csv
import json
import os
import re
from collections import defaultdict

import config as C

PVAR_RE = re.compile(r"\(p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\)")   # missense only
AA3 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E",
    "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F",
    "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
}
PLP = {"Pathogenic", "Likely Pathogenic",
       "Pathogenic, low penetrance", "Likely Pathogenic, low penetrance"}
DISTANCES = (0, 5, 10)
SUPP_MD = os.path.join(os.path.dirname(C.DATA_DIR), "..", "paper", "supp_pm1_reproducibility.md")
SUPP_MD = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                        "..", "paper", "supp_pm1_reproducibility.md"))
OUT_JSON = os.path.join(C.DATA_DIR, "pm1_reproducibility.json")


def has_pm1(codes):
    return any(t.strip().startswith("PM1") for t in (codes or "").split(","))


def load_variants():
    """One record per (variant identity, expert panel); track discordant re-curations."""
    seen = {}
    dup_discordant = []
    with open(C.EREPO_TSV) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            m = PVAR_RE.search(r.get("Variation", ""))
            if not m or m.group(1) not in AA3 or m.group(3) not in AA3:
                continue                                    # missense only (drop Ter/nonstandard)
            rec = {
                "g": (r.get("HGNC Gene Symbol") or "").strip(),
                "pos": int(m.group(2)),
                "ref": AA3.get(m.group(1), "?"), "alt": AA3.get(m.group(3), "?"),
                "panel": (r.get("Expert Panel") or "").strip(),
                "assert": (r.get("Assertion") or "").strip(),
                "met": has_pm1(r.get("Applied Evidence Codes (Met)")),
                "not_met": has_pm1(r.get("Applied Evidence Codes (Not Met)")),
            }
            vid = (r.get("ClinVar Variation Id") or "").strip() or r.get("Variation", "")
            key = (vid, rec["panel"])
            if key in seen:
                if seen[key]["met"] != rec["met"]:          # same variant+panel, PM1 flipped
                    dup_discordant.append((rec["g"], rec["panel"], rec["ref"], rec["pos"], rec["alt"]))
                continue
            seen[key] = rec
    return list(seen.values()), dup_discordant


def label(rec):
    return f"{rec['ref']}{rec['pos']}{rec['alt']}"


def find_pairs(variants, maxd, tier):
    """Discordant same-gene same-panel pairs. tier1 => both P/LP + explicit Not-Met."""
    by_gene = defaultdict(list)
    for v in variants:
        by_gene[v["g"]].append(v)
    seen, pairs = set(), []
    for g, rs in by_gene.items():
        rs = sorted(rs, key=lambda x: x["pos"])
        for i in range(len(rs)):
            for j in range(i + 1, len(rs)):
                a, b = rs[i], rs[j]
                if b["pos"] - a["pos"] > maxd:
                    break
                if a["panel"] != b["panel"] or not a["panel"]:
                    continue
                if tier == 1 and (a["assert"] not in PLP or b["assert"] not in PLP):
                    continue
                neg_b = b["not_met"] if tier == 1 else (not b["met"] and not b["not_met"])
                neg_a = a["not_met"] if tier == 1 else (not a["met"] and not a["not_met"])
                pos, neg = None, None
                if a["met"] and neg_b:
                    pos, neg = a, b
                elif b["met"] and neg_a:
                    pos, neg = b, a
                if pos is None:
                    continue
                key = (g, pos["panel"], min(a["pos"], b["pos"]), label(a), label(b))
                if key in seen:
                    continue
                seen.add(key)
                pairs.append({
                    "gene": g, "panel": pos["panel"], "dist": abs(a["pos"] - b["pos"]),
                    "a": {"aa": label(pos), "assertion": pos["assert"]},
                    "b": {"aa": label(neg), "assertion": neg["assert"]},
                    "tier": tier,
                })
    return pairs


def summarise(pairs):
    return {"pairs": len(pairs),
            "genes": len(set(p["gene"] for p in pairs)),
            "panels": len(set(p["panel"] for p in pairs))}


def main():
    variants, dup_discordant = load_variants()
    n_met = sum(v["met"] for v in variants)
    print(f"eRepo missense variants (deduped): {len(variants)} | PM1 met: {n_met} "
          f"| PM1 explicitly Not-Met: {sum(v['not_met'] for v in variants)}", flush=True)

    tier1 = {d: find_pairs(variants, d, 1) for d in DISTANCES}
    tier2 = {d: find_pairs(variants, d, 2) for d in DISTANCES}

    summary = {
        "n_missense": len(variants), "n_pm1_met": n_met,
        "tier1": {f"d{d}": summarise(tier1[d]) for d in DISTANCES},
        "tier2": {f"d{d}": summarise(tier2[d]) for d in DISTANCES},
        "same_variant_discordant": len(dup_discordant),
    }
    # headline examples: Tier-1, prefer same-residue then closest; a few well-known genes.
    ex = sorted(tier1[max(DISTANCES)], key=lambda p: (p["dist"], p["gene"]))
    examples = ex[:20]

    os.makedirs(C.DATA_DIR, exist_ok=True)
    with open(OUT_JSON, "w") as fh:
        json.dump({"summary": summary, "examples": examples}, fh, separators=(",", ":"))

    for d in DISTANCES:
        print(f"  Tier 1 (both P/LP, explicit Not-Met) Δ≤{d:2}: {summarise(tier1[d])}", flush=True)
    for d in DISTANCES:
        print(f"  Tier 2 (PM1 absent)                  Δ≤{d:2}: {summarise(tier2[d])}", flush=True)
    print(f"  Same variant, same panel, discordant PM1: {len(dup_discordant)}", flush=True)

    write_supp(summary, tier1, dup_discordant)
    print(f"Wrote {OUT_JSON}\nWrote {SUPP_MD}", flush=True)


def write_supp(summary, tier1, dup_discordant):
    lines = [
        "# Supplementary Table S1 — Inter-expert PM1 discordance in ClinGen eRepo",
        "",
        "*Generated by `scripts/analyze_pm1_reproducibility.py` from the ClinGen Evidence "
        "Repository export. Tier-1 pairs: same gene, same Expert Panel, both variants "
        "Pathogenic/Likely-Pathogenic, PM1 applied to variant A and **explicitly Not-Met** "
        "for variant B.*",
        "",
        f"Deduped eRepo missense variants: **{summary['n_missense']}** "
        f"(PM1 met: {summary['n_pm1_met']}). Tier-1 discordant pairs — "
        f"same residue: **{summary['tier1']['d0']['pairs']}** "
        f"({summary['tier1']['d0']['genes']} genes); "
        f"≤5 aa: **{summary['tier1']['d5']['pairs']}** "
        f"({summary['tier1']['d5']['genes']} genes); "
        f"≤10 aa: **{summary['tier1']['d10']['pairs']}** "
        f"({summary['tier1']['d10']['genes']} genes). "
        f"Same variant re-curated with discordant PM1: **{summary['same_variant_discordant']}**.",
        "",
        "| Gene | Expert Panel | PM1 applied | PM1 Not-Met | Δ aa |",
        "|---|---|---|---|---:|",
    ]
    for p in sorted(tier1[max(DISTANCES)], key=lambda x: (x["gene"], x["dist"])):
        lines.append(
            f"| {p['gene']} | {p['panel']} | {p['a']['aa']} ({p['a']['assertion']}) "
            f"| {p['b']['aa']} ({p['b']['assertion']}) | {p['dist']} |")
    with open(SUPP_MD, "w") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
