"""EXPLORATORY: do AlphaMissense islands land on *critical functional regions* (PM1's
"well-established functional domain"), or just generic conserved stretches?

We test, within island-bearing proteins, whether island residues are enriched — over
non-island residues of the same proteins — for three markers of functional criticality that
are each computed INDEPENDENTLY of AlphaMissense:

  1. UniProt functional domain   (curated structure/function)
  2. Cancer hotspot residue      (somatic recurrence, cancerhotspots.org; ±2 aa)
  3. gnomAD regional missense constraint  (population depletion; sub-regional o/e < 0.6)

Enrichment is a pooled odds ratio (island vs non-island residue), with a gene-cluster
bootstrap CI. Convergent enrichment across three orthogonal annotations argues islands
capture genuine functional criticality (PM1 semantics), not merely the predictor's own
conservation signal. Outputs paper/figures/fig_functional_enrichment.svg / .png.

Run:  python3 explore/island_functional_enrichment.py
"""
import glob
import gzip
import json
import os

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
GENES = os.path.join(os.path.expanduser("~"), "Dev", "am-islands-viewer", "public", "data", "genes")
FIG_BASE = os.path.normpath(os.path.join(HERE, "..", "paper", "figures", "fig_functional_enrichment"))
RMC_OE = 0.6
ANN = [("UniProt functional domain", "#2f7fa3"),
       ("Cancer hotspot (±2 aa)", "#b11949"),
       ("gnomAD constrained region", "#c98a2b")]


def masks(b):
    L = b["length"]
    isl = np.zeros(L, bool); dom = np.zeros(L, bool); rmc = np.zeros(L, bool); hot = np.zeros(L, bool)
    for i in (b.get("islands") or []):
        isl[max(0, i["s"] - 1):i["e"]] = True
    for d in (b.get("domains") or []):
        dom[max(0, d["s"] - 1):d["e"]] = True
    for r in (b.get("rmc") or []):
        if not r.get("whole") and r.get("oe", 1) < RMC_OE:
            rmc[max(0, r["s"] - 1):r["e"]] = True
    for h in (b.get("hotspots") or []):
        p = h["p"]; hot[max(0, p - 3):min(L, p + 2)] = True
    return isl, dom, rmc, hot


def main():
    files = glob.glob(os.path.join(GENES, "*.json.gz"))
    # per gene, per annotation: 2x2 [ [island&ann, island&~ann], [~island&ann, ~island&~ann] ]
    per_gene = []
    genes_with_island = 0
    for f in files:
        try:
            b = json.load(gzip.open(f, "rt"))
        except Exception:
            continue
        if not (b.get("islands")):
            continue
        genes_with_island += 1
        isl, dom, rmc, hot = masks(b)
        nout = (~isl)
        cells = {}
        for name, ann in (("dom", dom), ("hot", hot), ("rmc", rmc)):
            a = int((isl & ann).sum()); b_ = int((isl & ~ann).sum())
            c = int((nout & ann).sum()); d = int((nout & ~ann).sum())
            cells[name] = (a, b_, c, d)
        per_gene.append((b["gene"], cells))
    print(f"island-bearing genes: {genes_with_island}")

    def pooled_or(keys, sample):
        agg = {"dom": [0, 0, 0, 0], "hot": [0, 0, 0, 0], "rmc": [0, 0, 0, 0]}
        for idx in sample:
            _, cells = per_gene[idx]
            for k in agg:
                for j in range(4):
                    agg[k][j] += cells[k][j]
        out = {}
        for k in agg:
            a, b_, c, d = agg[k]
            out[k] = (a * d) / (b_ * c) if b_ and c else np.nan
        return out, agg

    rng = np.random.default_rng(20240715)
    n = len(per_gene)
    base, agg = pooled_or(["dom", "hot", "rmc"], range(n))
    boots = {"dom": [], "hot": [], "rmc": []}
    for _ in range(400):
        s = rng.integers(0, n, n)
        o, _ = pooled_or(None, s)
        for k in boots:
            boots[k].append(o[k])
    ci = {k: tuple(np.nanpercentile(boots[k], [2.5, 97.5])) for k in boots}

    # prevalence in island vs non-island residues
    prev = {}
    for k in ("dom", "hot", "rmc"):
        a, b_, c, d = agg[k]
        prev[k] = (a / (a + b_), c / (c + d))
    order = [("dom", ANN[0]), ("hot", ANN[1]), ("rmc", ANN[2])]
    print(f"{'annotation':28} {'OR':>6} {'95% CI':>15}  island%%  outside%%")
    for k, (name, _) in order:
        print(f"{name:28} {base[k]:6.2f} [{ci[k][0]:5.2f},{ci[k][1]:5.2f}] "
              f"{prev[k][0]*100:7.1f} {prev[k][1]*100:8.1f}")

    # ---------------------------------------------------------------- figure
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.linewidth": 0.7, "axes.edgecolor": "#3a4652", "svg.fonttype": "none"})
    fig, (axP, axO) = plt.subplots(1, 2, figsize=(11.2, 4.2), gridspec_kw={"width_ratios": [1.25, 1]})

    # Panel A: paired prevalence bars
    y = np.arange(len(order))[::-1]
    h = 0.34
    fmtpct = lambda p: (f"{p:.0f}%" if p >= 1 else f"{p:.2f}%")
    for i, (k, (name, col)) in enumerate(order):
        yi = y[i]
        axP.barh(yi + h / 2, prev[k][0] * 100, height=h, color=col, edgecolor="white", zorder=3)
        axP.barh(yi - h / 2, prev[k][1] * 100, height=h, color="#c7ccd4", edgecolor="white", zorder=3)
        axP.text(prev[k][0] * 100 + 0.6, yi + h / 2, fmtpct(prev[k][0] * 100), va="center", fontsize=8, color=col)
        axP.text(prev[k][1] * 100 + 0.6, yi - h / 2, fmtpct(prev[k][1] * 100), va="center", fontsize=8, color="#8999aa")
    axP.set_yticks(y); axP.set_yticklabels([nm for _, (nm, _) in order], fontsize=9)
    axP.set_xlabel("% of residues carrying the annotation")
    axP.set_title("A   Prevalence: island vs non-island residues", loc="left", fontsize=10.5,
                  fontweight="bold", color="#1c212c")
    axP.scatter([], [], marker="s", c="#2f7fa3", label="island residues")
    axP.scatter([], [], marker="s", c="#c7ccd4", label="non-island residues")
    axP.legend(loc="lower right", frameon=False, fontsize=8.3)

    # Panel B: OR forest (log)
    for i, (k, (name, col)) in enumerate(order):
        yi = y[i]
        lo, hi = ci[k]
        axO.plot([lo, hi], [yi, yi], color=col, lw=2, zorder=2)
        axO.scatter([base[k]], [yi], s=70, color=col, zorder=3, edgecolor="white", linewidth=1)
        axO.text(hi * 1.05, yi, f"{base[k]:.1f}× ({lo:.1f}–{hi:.1f})", va="center", fontsize=8.3, color="#1c212c")
    axO.axvline(1, color="#3a4652", lw=0.9)
    axO.set_xscale("log"); axO.set_yticks(y); axO.set_yticklabels([])
    axO.set_ylim(-0.6, len(order) - 0.4)
    axO.set_xlabel("Enrichment odds ratio (island vs non-island residue)")
    axO.set_title("B   Enrichment", loc="left", fontsize=10.5, fontweight="bold", color="#1c212c")

    fig.suptitle("AlphaMissense islands localise to independently-annotated functional regions",
                 x=0.5, y=1.0, fontsize=12, fontweight="bold", color="#1c212c")
    fig.text(0.5, -0.02, f"Within {genes_with_island:,} island-bearing proteins; enrichment vs the "
             f"same proteins' non-island residues (gene-cluster bootstrap 95% CI).",
             ha="center", fontsize=8, color="#8999aa")
    fig.tight_layout(rect=[0, 0.02, 1, 0.97])
    for ext in ("svg", "png"):
        fig.savefig(f"{FIG_BASE}.{ext}", dpi=300, bbox_inches="tight", facecolor="white")
    print("saved", FIG_BASE + ".svg / .png")


if __name__ == "__main__":
    main()
