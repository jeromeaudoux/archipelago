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

    # ---- figure: LR per stratum (log scale) with ACMG strength thresholds ----
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.linewidth": 0.7, "axes.edgecolor": "#3a4652",
                         "svg.fonttype": "none"})
    show = [r for r in rows if r[0] != "outside"]
    fig, ax = plt.subplots(figsize=(8.2, 4.2))
    xs = np.arange(len(show))
    lrs = np.array([r[1] for r in show])
    lo = np.array([r[2] for r in show]); hi = np.array([r[3] for r in show])
    cols = ["#8999aa" if r[0] == "island+benign" else "#b11949" for r in show]
    for name, y in THRESH.items():
        ax.axhline(y, ls="--", lw=0.8, color="#c7ccd4", zorder=0)
        ax.text(len(show) - 0.45, y * 1.02, name, fontsize=8, color="#8999aa", va="bottom", ha="right")
    ax.errorbar(xs, lrs, yerr=[lrs - lo, hi - lrs], fmt="none", ecolor="#c7ccd4", lw=1.4, zorder=2)
    ax.scatter(xs, lrs, s=46, c=cols, zorder=3, edgecolor="white", linewidth=0.7)
    for x, r in zip(xs, show):
        ax.annotate(f"{r[1]:.1f}", (x, r[3]), textcoords="offset points", xytext=(0, 5),
                    ha="center", fontsize=8, color="#1c212c")
    ax.axhline(1.0, color="#3a4652", lw=0.8)
    ax.set_yscale("log")
    ax.set_xticks(xs)
    ax.set_xticklabels([r[0].replace("island ", "").replace("+benign", "+ benign") for r in show],
                       rotation=20, ha="right", fontsize=8.5)
    ax.set_ylabel("Likelihood ratio for pathogenicity (log scale)")
    ax.set_xlabel("AM-island context (leave-one-out ClinVar corroboration)")
    ax.set_title("Graded PM1 calibration: in-island ClinVar corroboration lifts the likelihood "
                 "ratio\nfrom Supporting to Strong; benign-contradicted islands are downweighted",
                 fontsize=10, fontweight="bold", color="#1c212c")
    fig.tight_layout()
    for ext in ("svg", "png"):
        fig.savefig(f"{FIG_BASE}.{ext}", dpi=300, bbox_inches="tight", facecolor="white")
    print("\nsaved", FIG_BASE + ".svg / .png")


if __name__ == "__main__":
    main()
