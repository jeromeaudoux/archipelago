# Exploration — can a secondary rule rescue no-island PM1 false negatives?

*Exploratory only; nothing here is wired into the engine. Reproduce with
`python3 explore/rescue_no_island_pm1.py`. Base benchmark: TP 618, FP 266, FN 604
(P 0.699, R 0.506, F1 0.587).*

## Question

Island-based PM1 is high-precision but misses ~575 no-island false negatives (dominated by
sub-threshold ClinVar hotspots, §3.3). Is there a fallback rule — applied **only** where no
island overlaps — that flips true FN to TP without turning no-island true-negatives into
new FP?

## Setup

All eRepo missense residues the engine left PM1-off **because no island overlaps**:
4,510 residues → 3,724 mapped to a bundle (isoform-concordant). Of these **470 are FN**
(expert applied PM1) and **3,254 are TN** (expert did not) — a **~7:1 TN:FN base rate**.
That imbalance is the core obstacle.

## Signal separation (FN vs TN)

| Signal | FN | TN | Note |
|---|---|---|---|
| Residue AM | 0.914 | 0.664 | elevated in FN but TN also high |
| Windowed-mean AM (±5) | 0.712 | 0.589 | weak |
| gnomAD RMC o/e (sub-regional) | 0.737 | 0.819 | weak |
| Cancer hotspot ±2 aa | 35/470 (7.4 %) | 56/3254 (1.7 %) | best relative enrichment (~4×) |
| UniProt domain member | 219/470 (47 %) | 1067/3254 (33 %) | far too broad |

## Rescue-rule operating points

| Rule (no-island residues only) | rescued FN | new FP | P | R | F1 |
|---|--:|--:|--|--|--|
| raw AM residue > 0.5 | 407 | 1,987 | 0.313 | 0.839 | 0.456 |
| UniProt domain member | 219 | 1,067 | 0.386 | 0.685 | 0.494 |
| **[circular] ≥3 P/LP & 0 benign ±10** | 209 | 1,113 | 0.375 | 0.677 | 0.482 |
| gnomAD RMC o/e < 0.6 | 84 | 217 | 0.592 | 0.574 | 0.583 |
| just-outside-island ≤8 aa | 24 | 62 | 0.662 | 0.525 | 0.586 |
| cancer-hotspot ±2 aa | 35 | 56 | 0.670 | 0.534 | 0.594 |
| AM wmean>0.55 AND (hotspot OR RMC<0.5) | 91 | 161 | 0.624 | 0.580 | 0.601 |
| stacked safe (hotspot ∪ boundary ∪ AM&RMC) | 108 | 214 | 0.602 | 0.594 | 0.598 |

## Findings

1. **No clean rescue exists.** Every rule loses precision; the best F1 gain is trivial
   (0.587 → ~0.60). Since precision is the method's headline strength, trading 7–10 points
   of it for ~7 points of recall is a poor deal.

2. **The FN we most want are unreachable by orthogonal signals.** 345/470 mapped FN are
   sub-threshold ClinVar hotspots (≥3 P/LP nearby), but only 34 are cancer hotspots, 8 are
   near an island, and 70 are RMC-constrained. The remainder are "findable" only by ClinVar
   density itself — which floods FP (1,113 no-island TN also sit near ≥3 P/LP). Dense-ClinVar
   FN and dense-ClinVar TN occupy the same regions.

3. **Negative results reinforce the paper.** Raw-AM relaxation, domain membership, and the
   circular ClinVar-density rule all collapse precision to 0.31–0.39. This is empirical
   support for (a) why the island **contiguity + length gate** matters, and (b) the
   anti-circularity thesis — ClinVar density is a poor PM1 discriminator in this regime.

4. **The FP cost is an upper bound.** Because expert panels apply PM1 inconsistently to
   neighbours (§3.5), some TN→FP flips are arguably-correct calls against a noisy label; the
   true precision cost is likely smaller, but not cleanly quantifiable.

## If a rescue were ever adopted

The only defensible shape is a **conjunctive, orthogonal-evidence** rule on no-island
residues — sub-threshold-but-elevated AM **and** an independent corroborator (cancer-hotspot
recurrence or strong gnomAD regional constraint), explicitly **not** ClinVar density — and
applied only at **Supporting** strength via the existing graded framework, so a
low-confidence rescue does not over-commit. Even then it buys ≈ +9 recall points for ≈ −8
precision, so the recommendation is to leave it out: the island's specificity is the
differentiator and this is the wrong place to spend it.

## Possible paper use

A short Discussion paragraph ("why not a fallback rule for no-island residues") citing
finding 3 as further evidence for the island gate and against ClinVar-density PM1.
