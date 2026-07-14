"""Composed publication figure for the PM1 benchmark (Fig. 3).

Panel A — precision/recall of island-based PM1 vs Franklin and InterVar against expert
curation, with F1 iso-contours and the recalibrated operating point.
Panel B — the honest FP/FN divergence taxonomy (from public/data/pm1_benchmark.json).

Run:  python3 explore/paper_fig_benchmark.py   ->  paper/figures/fig_benchmark.svg / .png
"""
import json
import os
from collections import defaultdict

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
BENCH = os.path.join(HERE, "..", "public", "data", "pm1_benchmark.json")
FIG_BASE = os.path.normpath(os.path.join(HERE, "..", "paper", "figures", "fig_benchmark"))

# Published SeqOne benchmark (support = 1,270): (recall, precision, F1)
TOOLS = [
    ("AlphaMissense Islands", 0.50, 0.70, 0.59, "#b11949", "o"),
    ("Franklin", 0.86, 0.33, 0.47, "#0086e6", "s"),
    ("InterVar", 0.84, 0.19, 0.31, "#5b6673", "^"),
]
RECAL = (0.63, 0.84)  # illustrative recalibrated operating point (defensible FP / out-of-scope FN)

# FN / FP divergence categories: key -> (label, colour)
CAT = {
    "fn_subthreshold_hotspot": ("Sub-threshold ClinVar hotspot", "#b11949"),
    "fn_weak_am": ("Sparse evidence (non-AM)", "#e08a4e"),
    "fn_near_island": ("Just outside an island", "#edc765"),
    "fn_not_missense": ("Not missense (out of scope)", "#9aa3af"),
    "fp_clinvar_hotspot": ("ClinVar-supported hotspot (defensible)", "#2f7fa3"),
    "fp_am_only": ("AlphaMissense-only island", "#7fb2cc"),
    "fp_benign_conflict": ("Benign nearby (over-call)", "#b8d4e2"),
    "unmapped": ("Isoform-ambiguous", "#d3d8de"),
}
FN_ORDER = ["fn_subthreshold_hotspot", "fn_weak_am", "fn_near_island", "fn_not_missense", "unmapped"]
FP_ORDER = ["fp_clinvar_hotspot", "fp_am_only", "fp_benign_conflict", "unmapped"]


def main():
    data = json.load(open(BENCH))
    s = data["summary"]
    tally = defaultdict(lambda: defaultdict(int))
    for v in data["variants"]:
        if v["cat"] in ("fp", "fn"):
            tally[v["cat"]][v.get("dcat", "unmapped")] += 1
    fn, fp = tally["fn"], tally["fp"]
    NFN, NFP = sum(fn.values()), sum(fp.values())

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.linewidth": 0.7, "axes.edgecolor": "#3a4652", "svg.fonttype": "none"})
    fig = plt.figure(figsize=(12.4, 4.9))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.25], wspace=0.28)
    axA = fig.add_subplot(gs[0]); axB = fig.add_subplot(gs[1])

    # ---- Panel A: precision–recall with F1 iso-contours ----
    rr = np.linspace(0.001, 1, 400)
    for f1 in (0.3, 0.4, 0.5, 0.6, 0.7):
        with np.errstate(divide="ignore", invalid="ignore"):
            pp = (f1 * rr) / (2 * rr - f1)
        ok = (2 * rr > f1) & (pp <= 1.02)
        axA.plot(rr[ok], pp[ok], color="#d3d8de", lw=0.8, ls="--", zorder=0)
        rlab = 0.92                                    # label along each curve (well separated)
        plab = (f1 * rlab) / (2 * rlab - f1)
        axA.text(rlab, plab, f"F1={f1:g}", color="#aab0b8", fontsize=7,
                 ha="left", va="bottom", rotation=-30)
    for name, r, p, f1, col, mk in TOOLS:
        big = name.startswith("Alpha")
        axA.scatter([r], [p], s=150 if big else 90, c=col, marker=mk, zorder=5,
                    edgecolor="white", linewidth=1.2)
        dx, dy = (0.015, 0.03) if big else (0.0, -0.055)
        axA.annotate(f"{name}\nP={p:.2f} R={r:.2f} F1={f1:.2f}", (r, p),
                     xytext=(r + dx, p + dy), fontsize=8, color=col,
                     ha="left" if big else "center", va="bottom" if big else "top", fontweight="bold")
    # recalibrated point + arrow from AM islands
    axA.annotate("", xy=RECAL, xytext=(0.50, 0.70),
                 arrowprops=dict(arrowstyle="->", color="#b11949", lw=1.1, ls=(0, (3, 2))), zorder=4)
    axA.scatter([RECAL[0]], [RECAL[1]], s=95, facecolor="white", edgecolor="#b11949",
                linewidth=1.6, marker="o", zorder=5)
    axA.annotate("recalibrated\n(defensible FP / out-of-scope FN)", RECAL,
                 xytext=(RECAL[0] + 0.02, RECAL[1] + 0.02), fontsize=7.3, color="#b11949", va="bottom")
    axA.set_xlim(0, 1.02); axA.set_ylim(0, 1.06)
    axA.set_xlabel("Recall"); axA.set_ylabel("Precision")
    axA.set_title("A   PM1 vs expert curation", loc="left", fontsize=11, fontweight="bold", color="#1c212c")

    # ---- Panel B: divergence taxonomy stacked bars ----
    def draw_bar(y, counts, order, total, title):
        left = 0
        for k in order:
            n = counts.get(k, 0)
            if not n:
                continue
            lab, col = CAT[k]
            axB.barh(y, n, left=left, height=0.62, color=col, edgecolor="white", linewidth=1.2, zorder=3)
            if n / total > 0.06:
                axB.text(left + n / 2, y, f"{n}", ha="center", va="center", fontsize=8.5,
                         color="white" if k != "fp_benign_conflict" else "#33475b", fontweight="bold")
            left += n
        axB.text(-8, y, title, ha="right", va="center", fontsize=9.5, fontweight="bold", color="#1c212c")

    draw_bar(1, fn, FN_ORDER, NFN, f"False negatives\n(n={NFN})")
    draw_bar(0, fp, FP_ORDER, NFP, f"False positives\n(n={NFP})")
    axB.set_ylim(-0.6, 1.7); axB.set_xlim(0, max(NFN, NFP) * 1.02)
    axB.set_yticks([]); axB.spines["left"].set_visible(False)
    axB.tick_params(axis="y", length=0)
    axB.set_xlabel("variants")
    axB.set_title("B   Where the engine and ClinGen diverge", loc="left", fontsize=11,
                  fontweight="bold", color="#1c212c")
    # callouts
    axB.annotate(f"{fn['fn_subthreshold_hotspot']} real ClinVar hotspots below the AM threshold "
                 f"({fn['fn_subthreshold_hotspot']/NFN:.0%} of FN)",
                 xy=(fn["fn_subthreshold_hotspot"] / 2, 1.33), fontsize=7.6, color="#b11949", ha="center")
    axB.annotate(f"{fp['fp_clinvar_hotspot']} defensible, uncurated hotspots "
                 f"({fp['fp_clinvar_hotspot']/NFP:.0%} of FP)",
                 xy=(fp["fp_clinvar_hotspot"] / 2, -0.62), fontsize=7.6, color="#2f7fa3", ha="center")
    # legend
    handles = [plt.Rectangle((0, 0), 1, 1, color=CAT[k][1]) for k in
               ["fn_subthreshold_hotspot", "fn_weak_am", "fn_near_island", "fn_not_missense",
                "fp_clinvar_hotspot", "fp_am_only", "fp_benign_conflict", "unmapped"]]
    labels = [CAT[k][0] for k in
              ["fn_subthreshold_hotspot", "fn_weak_am", "fn_near_island", "fn_not_missense",
               "fp_clinvar_hotspot", "fp_am_only", "fp_benign_conflict", "unmapped"]]
    axB.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2,
               fontsize=7.3, frameon=False, handlelength=1.1, columnspacing=1.2)

    fig.suptitle("PM1 benchmark against ClinGen expert curation (eRepo)", x=0.5, y=1.02,
                 fontsize=12.5, fontweight="bold", color="#1c212c")
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    for ext in ("svg", "png"):
        fig.savefig(f"{FIG_BASE}.{ext}", dpi=300, bbox_inches="tight", facecolor="white")
    print(f"TP={s['tp']} FP={NFP} FN={NFN} | P={s['precision']} R={s['recall']} F1={s['f1']}")
    print("saved", FIG_BASE + ".svg / .png")


if __name__ == "__main__":
    main()
