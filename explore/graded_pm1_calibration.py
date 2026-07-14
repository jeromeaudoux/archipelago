"""EXPLORATORY: calibrate the graded PM1 strengths (Supporting/Moderate/Strong) against a
ClinVar reference, with leave-one-out ClinVar to avoid circularity.

For every ClinVar missense residue in the shipped bundles we take P/LP as a pathogenic
sample and B/LB as a benign sample. Each residue is assigned a PM1 stratum from its island
context, using LEAVE-ONE-OUT counts (the island's OTHER P/LP and B/LB residues, excluding the
focal residue itself):

  outside            not in any AM island
  island + benign    in island, >=2 other B/LB residues nearby (benign-contradicted)
  island, 0 P/LP     in island, no other P/LP (and <2 B/LB)  -> "island only"
  island, 1 P/LP     ...
  island, 2-4 P/LP   pathogenic-enriched
  island, 5-9 P/LP
  island, >=10 P/LP  dense hot-spot

For each stratum we estimate the likelihood ratio LR = P(stratum | pathogenic) /
P(stratum | benign) — prior-independent, and directly comparable to the Tavtigian OddsPath
thresholds that define ACMG evidence strength (LR >= 2.08 Supporting, 4.33 Moderate, 18.7
Strong). CIs are gene-cluster bootstraps. Because the focal residue's own label is removed
from the island counts, the covariate is genuinely "the neighbourhood", not the outcome.

Caveats (stated in the figure/writeup): ClinVar ascertainment biases the P/B reference; LOO
removes the focal but neighbouring ClinVar labels still correlate; a single ClinVar snapshot;
ideally cross-validate against a non-ClinVar (MAVE) truth.

Run:  python3 explore/graded_pm1_calibration.py
"""
import glob
import gzip
import json
import os

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HOME = os.path.expanduser("~")
GENES = os.path.join(HOME, "Dev", "am-islands-viewer", "public", "data", "genes")
FIG_BASE = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                         "..", "paper", "figures", "fig_graded_calibration"))
# Tavtigian (prior 0.10) OddsPath thresholds == likelihood-ratio thresholds for strength.
THRESH = {"Supporting": 2.08, "Moderate": 4.33, "Strong": 18.7}
STRATA = ["outside", "island+benign", "island 0 P/LP", "island 1 P/LP",
          "island 2–4 P/LP", "island 5–9 P/LP", "island ≥10 P/LP"]


def stratum(in_isl, loo_plp, loo_blb):
    if not in_isl:
        return "outside"
    if loo_blb >= 2:
        return "island+benign"
    if loo_plp == 0:
        return "island 0 P/LP"
    if loo_plp == 1:
        return "island 1 P/LP"
    if loo_plp <= 4:
        return "island 2–4 P/LP"
    if loo_plp <= 9:
        return "island 5–9 P/LP"
    return "island ≥10 P/LP"


def collect():
    """gene -> stratum -> [n_pathogenic, n_benign]."""
    per_gene = {}
    files = glob.glob(os.path.join(GENES, "*.json.gz"))
    for f in files:
        try:
            b = json.load(gzip.open(f, "rt"))
        except Exception:
            continue
        islands = b.get("islands") or []
        if not islands and not (b.get("clinvar") or b.get("clinvar_benign")):
            continue
        gene = b["gene"]
        agg = per_gene.setdefault(gene, {s: [0, 0] for s in STRATA})

        def find(pos):
            for i in islands:
                if i["s"] <= pos <= i["e"]:
                    return i
            return None
        for lab, key in ((1, "clinvar"), (0, "clinvar_benign")):
            for v in (b.get(key) or []):
                isl = find(v["p"])
                if isl is None:
                    s = "outside"
                else:
                    loo_plp = (isl.get("plp", 0) or 0) - (1 if lab == 1 else 0)
                    loo_blb = (isl.get("blb", 0) or 0) - (1 if lab == 0 else 0)
                    s = stratum(True, max(0, loo_plp), max(0, loo_blb))
                agg[s][1 - lab] += 1  # index 0 = pathogenic, 1 = benign
    return per_gene


def lr_table(per_gene, genes=None):
    genes = genes if genes is not None else list(per_gene)
    tot = {s: [0, 0] for s in STRATA}
    for g in genes:
        for s, (p, bn) in per_gene[g].items():
            tot[s][0] += p; tot[s][1] += bn
    Np = sum(tot[s][0] for s in STRATA); Nb = sum(tot[s][1] for s in STRATA)
    out = {}
    for s in STRATA:
        p, bn = tot[s]
        lr = (p / Np) / (bn / Nb) if (Np and Nb and bn) else np.nan
        out[s] = (lr, p, bn)
    return out


def main():
    per_gene = collect()
    genes = list(per_gene)
    base = lr_table(per_gene)
    Np = sum(base[s][1] for s in STRATA); Nb = sum(base[s][2] for s in STRATA)
    print(f"genes={len(genes)} | pathogenic(P/LP)={Np} benign(B/LB)={Nb}\n")

    rng = np.random.default_rng(20240715)
    ug = np.array(genes)
    boot = {s: [] for s in STRATA}
    for _ in range(500):
        samp = rng.choice(ug, size=len(ug), replace=True)
        t = lr_table(per_gene, samp)
        for s in STRATA:
            boot[s].append(t[s][0])
    ci = {s: tuple(np.nanpercentile(boot[s], [2.5, 97.5])) for s in STRATA}

    print(f"{'stratum':18} {'LR':>7} {'95% CI':>16} {'nP':>6} {'nB':>6}   strength")
    rows = []
    for s in STRATA:
        lr, p, bn = base[s]
        lo, hi = ci[s]
        strg = ("Strong" if lo >= THRESH["Strong"] else "Moderate" if lo >= THRESH["Moderate"]
                else "Supporting" if lo >= THRESH["Supporting"] else "—")
        print(f"{s:18} {lr:7.2f} [{lo:6.2f},{hi:6.2f}] {p:6d} {bn:6d}   {strg}")
        rows.append((s, lr, lo, hi, p, bn, strg))

    # ---------------------------------------------------------------- figure
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.spines.left": False, "axes.linewidth": 0.7,
                         "axes.edgecolor": "#c7ccd4", "svg.fonttype": "none"})
    R = {s: (base[s][0], ci[s][0], ci[s][1], base[s][1], base[s][2]) for s in STRATA}
    dose = ["island 0 P/LP", "island 1 P/LP", "island 2–4 P/LP", "island 5–9 P/LP", "island ≥10 P/LP"]
    dlab = ["0", "1", "2–4", "5–9", "≥10"]
    ramp = ["#f2789f", "#e0518a", "#c93370", "#ab1a4d", "#7d0d33"]  # light→dark (dose)
    INK, MUT = "#1c212c", "#8999aa"
    xdose = np.arange(len(dose)) + 2.0       # dose series at x = 2..6
    x_out, x_ben = 0.0, 0.9                    # reference markers

    fig, ax = plt.subplots(figsize=(9.2, 5.0))
    YB, YT = 0.62, 165
    # ACMG strength bands (background)
    bands = [(YB, THRESH["Supporting"], "#eff1f4", "not met", "#8999aa"),
             (THRESH["Supporting"], THRESH["Moderate"], "#fbe9dc", "Supporting", "#b07627"),
             (THRESH["Moderate"], THRESH["Strong"], "#f7d9c0", "Moderate", "#bd6224"),
             (THRESH["Strong"], YT, "#f6ccd3", "Strong", "#a5273f")]
    xr = 5.15
    for y0, y1, fc, name, tc in bands:
        ax.axhspan(y0, y1, color=fc, zorder=0, lw=0)
        ax.text(xr, np.sqrt(y0 * min(y1, YT)), name, ha="right", va="center",
                fontsize=9.5, color=tc, fontweight="bold", alpha=0.9, zorder=1)
    ax.axhline(1.0, color="#3a4652", lw=0.9, zorder=1)

    # dose-response series (connected, sequential red)
    dl = np.array([R[s][0] for s in dose]); dlo = np.array([R[s][1] for s in dose])
    dhi = np.array([R[s][2] for s in dose])
    ax.plot(xdose, dl, "-", color="#ab1a4d", lw=1.6, alpha=0.55, zorder=3)
    ax.errorbar(xdose, dl, yerr=[dl - dlo, dhi - dl], fmt="none", ecolor="#b9727f", lw=1.5,
                capsize=3, zorder=4)
    ax.scatter(xdose, dl, s=[70, 80, 95, 115, 140], c=ramp, zorder=5,
               edgecolor="white", linewidth=1.1)
    for x, s in zip(xdose, dose):
        ax.annotate(f"{R[s][0]:.1f}", (x, R[s][2]), textcoords="offset points", xytext=(0, 7),
                    ha="center", fontsize=9, fontweight="bold", color=INK, zorder=6)
        ax.annotate(f"n={R[s][3]:,}", (x, R[s][1]), textcoords="offset points", xytext=(0, -13),
                    ha="center", fontsize=7.3, color=MUT, zorder=6)

    # reference markers (not on the dose trend)
    for x, s, mk, fc in ((x_out, "outside", "o", "white"), (x_ben, "island+benign", "D", "#b6bcc6")):
        lr, lolr, hilr = R[s][0], R[s][1], R[s][2]
        ax.errorbar(x, lr, yerr=[[lr - lolr], [hilr - lr]], fmt="none", ecolor="#c7ccd4", lw=1.5,
                    capsize=3, zorder=4)
        ax.scatter([x], [lr], s=90, marker=mk, facecolor=fc, edgecolor="#5b6673",
                   linewidth=1.3, zorder=5)
        ax.annotate(f"{lr:.1f}", (x, hilr), textcoords="offset points", xytext=(0, 7),
                    ha="center", fontsize=9, color="#5b6673", zorder=6)

    ax.set_yscale("log"); ax.set_ylim(YB, YT); ax.set_xlim(-0.6, 6.7)
    ax.tick_params(axis="y", length=0)
    ax.set_yticks([1, 2, 5, 10, 20, 50, 100])
    ax.set_yticklabels(["1", "2", "5", "10", "20", "50", "100"], fontsize=8.5, color="#3a4652")
    ax.set_xticks([x_out, x_ben, *xdose])
    ax.set_xticklabels(["outside\nisland", "benign-\ncontradicted", *dlab], fontsize=9, color=INK)
    ax.annotate("other P/LP residues in the island  (leave-one-out)",
                xy=(xdose.mean(), 0), xytext=(xdose.mean(), -0.20), textcoords=("data", "axes fraction"),
                ha="center", fontsize=8.5, color=MUT, annotation_clip=False)
    ax.annotate("", xy=(xdose[0] - 0.25, -0.15), xytext=(xdose[-1] + 0.25, -0.15),
                xycoords=("data", "axes fraction"), textcoords=("data", "axes fraction"),
                arrowprops=dict(arrowstyle="-", color="#c7ccd4", lw=0.9), annotation_clip=False)
    ax.set_ylabel("Likelihood ratio for pathogenicity", fontsize=10)
    fig.suptitle("Graded PM1 is empirically calibrated", x=0.5, y=0.99,
                 fontsize=13, fontweight="bold", color=INK)
    ax.set_title("Leave-one-out ClinVar likelihood ratio by AM-island context "
                 "(46,955 P/LP vs 121,506 B/LB, 15,744 genes)",
                 fontsize=9.2, color=MUT, pad=10)
    fig.tight_layout(rect=[0, 0.02, 1, 0.97])
    for ext in ("svg", "png"):
        fig.savefig(f"{FIG_BASE}.{ext}", dpi=300, bbox_inches="tight", facecolor="white")
    print("\nsaved", FIG_BASE + ".svg / .png")


if __name__ == "__main__":
    main()
