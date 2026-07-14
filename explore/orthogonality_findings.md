# Orthogonality — does an island add signal beyond the per-variant predictor (PP3)?

*Reproduce: `python3 explore/orthogonality_pp3_island.py` (needs the SeqOne scored ClinVar
JSON + the island BED). Figure: `paper/figures/fig_orthogonality.png`.*

## The question (falsifiable by design)

At a **fixed per-variant predictor score**, does being **inside an AlphaMissense island** still
change the odds of pathogenicity? If yes, the regional signal is non-redundant and co-applying
PM1 (island) with PP3 (the per-variant score) is not double-counting. If no, the island is just
a smoothed restatement of the per-residue score and the graded-PM1 story collapses. We run it
against three predictors and report all three — the AlphaMissense conditioning can falsify the
claim, which is what makes it credible.

## Data & method

ClinVar-labelled missense variants with populated AlphaMissense / REVEL / CADD from the SeqOne
scored set: **N = 3,133** (1,244 P/LP, 1,889 B/LB; 359 in-island, 2,774 out). Island membership
is **binary genomic overlap** with the island BED (transcript-agnostic, as the engine applies
it). Primary result = stratified pathogenic rate for in- vs out-of-island within predictor
deciles + a Cochran–Mantel–Haenszel common OR (conditions on the score) with a **gene-cluster
bootstrap CI**. Support = a **gene-cluster-robust** logistic regression
`pathogenic ~ z(score) + in_island`. Because `in_island` is binary membership (not
`island_mean_score`), the focal residue never enters the covariate → **VIF ≈ 1**, sidestepping
the leave-one-out collinearity trap.

## Result

| Conditioned on | CMH OR (in-island) | 95% gene-cluster CI | logit β | p | VIF |
|---|--:|---|--:|--:|--:|
| **REVEL** | **4.97** | 3.05 – 8.55 | 1.63 | 7×10⁻⁹ | 1.12 |
| **CADD (phred)** | **5.04** | 3.27 – 7.99 | 1.65 | 1×10⁻¹⁴ | 1.08 |
| AlphaMissense (per-variant) | 1.33 | 0.77 – 2.34 | 0.48 | 0.09 (n.s.) | 1.22 |

## Interpretation

1. **Against REVEL and CADD — strongly non-redundant.** At the same REVEL/CADD score, an
   in-island variant has ~**5× the odds** of being pathogenic. Because REVEL/CADD are the PP3
   sources actually deployed, this **licenses co-applying PM1 (island) with PP3 — it is not
   double-counting.** This is the result that defuses the SVI double-counting objection.

2. **Against AlphaMissense itself — redundant (as expected).** Conditioned on the variant's own
   AM score, island membership adds no significant signal (OR 1.3, p = 0.09; the two lines in the
   figure overlap). The island *is* essentially a smoothed restatement of the per-residue AM
   score — which is exactly why it doubles as the islands-vs-plain-threshold ablation. Low VIF
   (1.22) shows this null is genuine, not a collinearity artifact.

3. **The co-application rule that follows.** PM1-from-AM-islands and PP3 may be co-applied **iff
   PP3 is not itself AlphaMissense.** With a REVEL/CADD-based PP3 (SeqOne's deployment) the two
   are orthogonal and stack legitimately; pairing AM-island PM1 with an **AM-based PP3** would
   double-count and must be avoided (or capped). State this explicitly in the engine and paper.

## Caveats

- N = 3,133 is the ClinVar subset with all three predictors populated; a larger set (e.g.
  pooling `erepo_variants_updated2`) would tighten CIs. In-island variants are sparse at low
  predictor deciles (islands are high-AM regions), hence the wide in-island band there.
- The truth label is ClinVar clinical significance. The *conditional* comparison (in- vs
  out-of-island at fixed score) is not circular with respect to islands, since islands are
  AlphaMissense-derived and independent of ClinVar; the marginal base rates are not the claim.

## Why we do NOT add an "in-island × ClinVar overlap" stratum here

Tempting, but circular for *this* test: the outcome is ClinVar clinical significance, so a
covariate built from nearby ClinVar pathogenic variants would use ClinVar to predict ClinVar,
inflating the OR and re-importing the exact circularity islands were designed to avoid. Fig. 5
must stay a clean AlphaMissense-only-vs-PP3 comparison. "In-island **+ ClinVar corroboration**"
belongs to the *graded-PM1 calibration* (§4 / C3), and even there it is only defensible with
**leave-one-out ClinVar** (exclude the focal variant's own ClinVar status from the island's
density) and ideally a non-ClinVar truth label (e.g. MAVE). That is a separate analysis, not
part of the orthogonality claim.

## Output

`paper/figures/fig_orthogonality.svg` (vector, editable text) and `.png` (300 dpi). Caption
note: VIF ≈ 1.1 for all three conditionings (binary membership ≠ `island_mean_score`), so the
regression is stable and the AlphaMissense null is not a collinearity artifact; decile points
with < 10 variants in a line are suppressed.

## Paper use

Fills Results §3.4 (orthogonality) and answers open item C2. Headline for the Discussion:
islands add PP3-relevant regional signal over REVEL/CADD (co-application safe), but not over
AlphaMissense itself (so PP3 must not be AM when PM1 comes from AM islands).
