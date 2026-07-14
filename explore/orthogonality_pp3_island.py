"""Does island membership add pathogenicity signal BEYOND the per-variant predictor (PP3)?

Falsifiable test (§3.4 / C2): at a fixed per-variant score, does being inside an AlphaMissense
island still change the odds of pathogenicity?
  - If YES against REVEL/CADD (the PP3 sources actually deployed), the regional signal is
    non-redundant → co-applying PM1 (island) with PP3 (REVEL/CADD) is not double-counting.
  - If NO against AlphaMissense itself, the island is a smoothed restatement of the per-residue
    AM score → PM1-from-AM-islands would double-count with an AM-based PP3.

We report BOTH conditionings. The primary, collinearity-robust result is a stratified
comparison (pathogenic rate for in-island vs out-of-island within predictor deciles) plus a
Cochran–Mantel–Haenszel common odds ratio with a gene-cluster bootstrap CI. A cluster-robust
(by gene) logistic regression `pathogenic ~ z(score) + in_island` and the VIF are reported as
support. `in_island` is BINARY genomic membership (not island_mean_score), so the focal residue
does not enter the covariate — VIF stays ~1, sidestepping the leave-one-out collinearity trap.

Data (SeqOne, not in this repo): clinvar_acmg_small_updated.json.gz (ClinVar-labelled, with
AlphaMissense/REVEL/CADD) and data/am_islands_hg38.bed.gz. Outputs paper/figures/fig_orthogonality.png.

Run:  python3 explore/orthogonality_pp3_island.py
"""
import bisect
import gzip
import json
import os
import pickle

import matplotlib
import numpy as np
from scipy import stats

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HOME = os.path.expanduser("~")
DATA = os.path.join(HOME, "Dev", "acmgscore-v2", "clinvar_acmg_small_updated.json.gz")
BED = os.path.join(HOME, "Dev", "acmgscore-v2", "data", "am_islands_hg38.bed.gz")
CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".orth_cache.pkl")
FIG_BASE = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                         "..", "paper", "figures", "fig_orthogonality"))
MIN_N = 10  # suppress decile points with fewer than this many variants in a line
PATHO = {"Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"}
BENIGN = {"Benign", "Likely_benign", "Benign/Likely_benign"}


def load_islands():
    iv = {}
    with gzip.open(BED, "rt") as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            c = ln.split("\t")
            iv.setdefault(c[0], []).append((int(c[1]), int(c[2])))
    for ch in iv:
        iv[ch].sort()
    starts = {ch: [s for s, _ in v] for ch, v in iv.items()}
    def in_island(chrom, pos):
        ch = str(chrom) if str(chrom).startswith("chr") else "chr" + str(chrom)
        v = iv.get(ch)
        if not v:
            return False
        i = bisect.bisect_right(starts[ch], pos - 1) - 1
        return i >= 0 and v[i][0] <= pos - 1 < v[i][1]
    return in_island


def num(x):
    if isinstance(x, list):
        xs = [v for v in x if isinstance(v, (int, float))]
        return max(xs) if xs else None
    return x if isinstance(x, (int, float)) else None


def build_rows():
    if os.path.exists(CACHE):
        return pickle.load(open(CACHE, "rb"))
    in_island = load_islands()
    rows = []
    with gzip.open(DATA, "rt") as fh:
        for ln in fh:
            if not ln.strip():
                continue
            e = json.loads(ln)["element"]
            tags = e.get("extdb", {}).get("clinvar_tags") or []
            tag = tags[0] if len(tags) == 1 else None
            lab = 1 if tag in PATHO else 0 if tag in BENIGN else None
            if lab is None:
                continue
            t = (e.get("annotation") or {}).get("transcripts") or []
            tr = next((x for x in t if x.get("fsc_sc_am") is not None), None)  # missense proxy
            pos = e.get("estart")
            if tr is None or pos is None:
                continue
            g = (tr.get("ensembl_geneid") or "").split("&")[0].upper() or f"{e.get('echr')}:{pos}"
            rows.append((lab, num(tr.get("fsc_sc_am")), num(tr.get("revel_sc")),
                         num(tr.get("fsc_sc_cadd_phred")),
                         1 if in_island(e.get("echr"), int(pos)) else 0, g))
    pickle.dump(rows, open(CACHE, "wb"))
    return rows


def wilson(k, n):
    if n == 0:
        return (np.nan, np.nan)
    z, p, d = 1.96, k / n, 1 + 1.96 ** 2 / n
    c = (p + 1.96 ** 2 / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + 1.96 ** 2 / (4 * n * n)) / d
    return (c - h, c + h)


def cmh_or(y, isl, strat):
    num_, den = 0.0, 0.0
    for s in np.unique(strat):
        m = strat == s
        a = ((isl == 1) & (y == 1) & m).sum(); b = ((isl == 1) & (y == 0) & m).sum()
        c = ((isl == 0) & (y == 1) & m).sum(); dd = ((isl == 0) & (y == 0) & m).sum()
        n = a + b + c + dd
        if n:
            num_ += a * dd / n; den += b * c / n
    return num_ / den if den > 0 else np.nan


def deciles(x, q=10):
    edges = np.unique(np.nanquantile(x, np.linspace(0, 1, q + 1)))
    return np.clip(np.digitize(x, edges[1:-1]), 0, len(edges) - 2), edges


def irls_logit(X, y, ridge=1e-6, it=100):
    b = np.zeros(X.shape[1])
    p = np.full(len(y), 0.5)
    for _ in range(it):
        p = 1 / (1 + np.exp(-(X @ b)))
        W = np.clip(p * (1 - p), 1e-6, None)
        step = np.linalg.solve((X.T * W) @ X + ridge * np.eye(X.shape[1]), X.T @ (y - p) - ridge * b)
        b = b + step
        if np.max(np.abs(step)) < 1e-8:
            break
    return b, p


def cluster_robust_se(X, y, p, groups):
    W = np.clip(p * (1 - p), 1e-6, None)
    bread = np.linalg.inv((X.T * W) @ X)
    u = X * ((y - p)[:, None])
    meat = sum(np.outer(u[groups == g].sum(0), u[groups == g].sum(0)) for g in np.unique(groups))
    return np.sqrt(np.diag(bread @ meat @ bread))


def main():
    rows = build_rows()
    lab = np.array([r[0] for r in rows], float)
    isl = np.array([r[4] for r in rows], float)
    gene = np.array([r[5] for r in rows])
    preds = {"AlphaMissense (per-variant)": np.array([np.nan if r[1] is None else r[1] for r in rows]),
             "REVEL": np.array([np.nan if r[2] is None else r[2] for r in rows]),
             "CADD (phred)": np.array([np.nan if r[3] is None else r[3] for r in rows])}
    print(f"N={len(rows)} | patho={int(lab.sum())} benign={int(len(lab)-lab.sum())} "
          f"| in_island={int(isl.sum())} out={int(len(isl)-isl.sum())}")

    rng = np.random.default_rng(20240714)
    IN, OUT, INK, MUT = "#b11949", "#0086e6", "#1c212c", "#8999aa"
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9, "axes.titlesize": 10.5,
        "axes.labelsize": 9.5, "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": 0.7, "axes.edgecolor": "#3a4652",
        "xtick.labelsize": 8.5, "ytick.labelsize": 8.5, "xtick.color": "#3a4652",
        "ytick.color": "#3a4652", "legend.fontsize": 8.5, "svg.fonttype": "none",
    })

    def fmt_p(p):
        if p >= 1e-3:
            return f"p = {p:.3f}"
        e = int(np.floor(np.log10(p))); mant = p / 10 ** e
        sup = str(e).translate(str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))
        return f"p = {mant:.0f}×10{sup}"

    fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.7), sharey=True)
    for pi, (ax, (name, x)) in enumerate(zip(axes, preds.items())):
        m = ~np.isnan(x)
        xx, yy, ii, gg = x[m], lab[m], isl[m], gene[m]
        dec, edges = deciles(xx)
        cx, lines = [], {1: ([], []), 0: ([], [])}
        for k in range(len(edges) - 1):
            s = dec == k
            cx.append(np.median(xx[s]) if s.any() else np.nan)
            for grp in (1, 0):
                sg = s & (ii == grp); n = int(sg.sum()); kk = int(yy[sg].sum())
                lines[grp][0].append(kk / n if n >= MIN_N else np.nan)
                lines[grp][1].append(wilson(kk, n) if n >= MIN_N else (np.nan, np.nan))
        cx = np.array(cx)
        ax.grid(axis="y", color="#eef1f5", lw=0.8, zorder=0)
        for grp, col, ls, lbl in ((0, OUT, "--", "outside island"), (1, IN, "-", "in island")):
            p = np.array(lines[grp][0], float); ok = ~np.isnan(p)
            lo = np.array([a for a, _ in lines[grp][1]]); hi = np.array([b for _, b in lines[grp][1]])
            ax.fill_between(cx[ok], lo[ok], hi[ok], color=col, alpha=0.12, lw=0, zorder=1)
            ax.plot(cx[ok], p[ok], ls, color=col, lw=1.8, marker="o", ms=4.5, mec="white",
                    mew=0.6, label=lbl, zorder=3)
        orr = cmh_or(yy, ii, dec)
        ug = np.unique(gg)
        boots = [cmh_or(yy[idx], ii[idx], dec[idx]) for idx in
                 (np.concatenate([np.where(gg == s)[0] for s in rng.choice(ug, len(ug), replace=True)])
                  for _ in range(600))]
        lo_ci, hi_ci = np.nanpercentile(boots, [2.5, 97.5])
        z = (xx - np.nanmean(xx)) / np.nanstd(xx)
        X = np.column_stack([np.ones_like(z), z, ii])
        b, p = irls_logit(X, yy)
        se = cluster_robust_se(X, yy, p, gg)[2]
        beta, pval = b[2], 2 * stats.norm.sf(abs(b[2] / se))
        vif = 1 / (1 - np.corrcoef(z, ii)[0, 1] ** 2)
        ax.set_title(name, fontweight="bold", color=INK, pad=8)
        ax.set_xlabel(name)
        ax.set_ylim(-0.02, 1.04)
        ax.text(-0.04, 1.10, "ABC"[pi], transform=ax.transAxes, fontsize=13,
                fontweight="bold", va="top", ha="right", color=INK)
        sig = pval < 0.05
        ax.text(0.04, 0.96,
                f"OR = {orr:.1f}  ({lo_ci:.1f}–{hi_ci:.1f})\n{fmt_p(pval)}"
                + ("" if sig else "\n(n.s.)"),
                transform=ax.transAxes, va="top", fontsize=8.5,
                color=IN if sig else MUT,
                bbox=dict(boxstyle="round,pad=0.35", fc="white",
                          ec=IN if sig else "#dbe0e6", lw=0.8))
        print(f"  {name:26} CMH_OR={orr:.2f} [{lo_ci:.2f}-{hi_ci:.2f}]  beta={beta:.3f} p={pval:.2e} "
              f"VIF={vif:.2f} N={int(m.sum())}")
    axes[0].set_ylabel("P(pathogenic), within predictor decile")
    axes[0].legend(loc="lower right", frameon=False, handlelength=1.8)
    fig.suptitle("At a fixed per-variant score, island membership still raises pathogenicity odds",
                 fontsize=11.5, fontweight="bold", color=INK, y=1.02)
    fig.tight_layout()
    for ext in ("svg", "png"):
        fig.savefig(f"{FIG_BASE}.{ext}", dpi=300, bbox_inches="tight",
                    facecolor="white", transparent=False)
    print("saved", FIG_BASE + ".svg / .png")


if __name__ == "__main__":
    main()
