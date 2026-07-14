# Graded-PM1 calibration — likelihood ratios by island context (leave-one-out ClinVar)

*Reproduce: `python3 explore/graded_pm1_calibration.py`. Figure:
`paper/figures/fig_graded_calibration.svg` / `.png`.*

## Question

Does the proposed PM1 grading (Supporting → Moderate → Strong from island membership +
in-island ClinVar corroboration) hold up empirically, and what evidence strength does each
level actually earn?

## Method

Reference set from the shipped bundles: every ClinVar missense residue is a sample —
**P/LP = pathogenic (n = 46,955), B/LB = benign (n = 121,506)** across 15,744 genes. Each
residue is assigned a PM1 stratum from its island context using **leave-one-out ClinVar**
(the island's P/LP and B/LB counts with the focal residue removed), so the covariate is the
neighbourhood, never the outcome. For each stratum we estimate the **likelihood ratio**
LR = P(stratum | pathogenic) / P(stratum | benign) — prior-independent and directly comparable
to the Tavtigian OddsPath thresholds that define ACMG strength (LR ≥ 2.08 Supporting, ≥ 4.33
Moderate, ≥ 18.7 Strong). 95% CIs are gene-cluster bootstraps.

## Result

| Island context (LOO) | LR | 95% CI | nP | nB | Strength (by CI lower bound) |
|---|--:|---|--:|--:|---|
| outside island | 0.77 | 0.74–0.79 | 35,281 | 118,824 | — (not met) |
| island, benign-contradicted (≥2 B/LB) | 4.81 | 3.80–6.33 | 1,838 | 988 | Supporting |
| **island, 0 other P/LP** | **2.50** | 2.20–2.88 | 893 | 923 | **Supporting** |
| island, 1 other P/LP | 7.91 | 6.65–9.90 | 724 | 237 | Moderate |
| island, 2–4 P/LP | 15.64 | 13.37–18.49 | 1,656 | 274 | Moderate |
| island, 5–9 P/LP | 36.58 | 29.42–46.98 | 1,894 | 134 | Strong |
| island, ≥10 P/LP | 95.89 | 76.29–119.14 | 4,669 | 126 | Strong |

## Interpretation

- **An AM island with no ClinVar corroboration already reaches Supporting** (LR 2.5,
  CI 2.2–2.9 ≥ 2.08). This is the empirical basis for triggering PM1 at Supporting by default
  on the ~6,150 island-bearing genes that lack a VCEP and ClinVar — the paper's headline claim.
- **In-island ClinVar corroboration produces a clean dose-response**: 1 → 2–4 → 5–9 → ≥10 other
  P/LP residues lift the LR 7.9 → 15.6 → 36.6 → 95.9, crossing Moderate and then Strong. The
  proposed Supporting/Moderate/Strong ladder is therefore empirically supported, not heuristic.
- **Benign variation downweights**: islands with ≥2 other B/LB residues fall back to LR ~4.8
  (Supporting), well below the clean P-enriched levels — matching the "downweight when
  benign-contradicted" rule.
- **Outside islands, LR < 1** (0.77): PM1-not-met is the correct default.

The strengths map almost exactly onto the framework in the paper's grading table, moving it
from "proposed" to "ClinVar-calibrated (leave-one-out)".

## Caveats

- ClinVar ascertainment biases the P/B reference (this is the same reference class ClinGen
  PP3/BP4 calibration uses, but it is not a random sample of all variants).
- LOO removes the focal residue, but neighbouring ClinVar labels still correlate within a
  region — some of the corroboration LR reflects that genuine regional clustering (which is
  what PM1 encodes) rather than an independent axis.
- A single ClinVar snapshot (v107); strengths should be re-derived on update and ideally
  cross-validated against a non-ClinVar truth (e.g. MAVE functional scores) to fully break the
  residual circularity.

## Paper use

Fills §4 / open item C3 (Fig. 6). Note that this calibrates PM1 *strength*; §3.4/Fig. 5
already established that the criterion is non-redundant with the PP3 predictor it co-applies
with (REVEL/CADD).
