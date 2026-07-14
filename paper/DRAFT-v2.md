# AlphaMissense Islands: a genome-wide resource of regional missense hotspots for non-circular, reproducible ACMG PM1

*Advanced working draft (v2) — SeqOne / Archipelago. Live tool:
[archipelago2.vercel.app](https://archipelago2.vercel.app). Shareable Google Doc:
[docs.google.com/…/AlphaMissense-Islands-v2](https://docs.google.com/document/d/1uKvgLb12fCJefmFWL2pxy2e06LT_9XmgJw3y_EN68xs/edit).
Author list TBD (J. Audoux, J. B. [surname], et al.). Citation keys `[CITE:key]` resolve to
the Reference list; `[verify]` marks IDs to confirm at submission. Figures reference
`paper/figures/`. `DRAFT.md` is the previous version.*

---

## Target format & venue

Short **Application Note / Brief Communication** (≈2,000–2,800 words, 3 main figures,
2 tables). Path: **bioRxiv/medRxiv preprint first** (priority vs. the gnomAD
regional-constraint→PM1 preprint), then **Bioinformatics Advances** or **Bioinformatics**
(tooling audience), **NAR Web Server** (resource/reach), or **HGG Advances / GIM Open**
(clinical audience).

---

## Abstract

The ACMG/AMP criterion **PM1** — a variant located in a mutational hotspot or a
critical, well-established functional domain — is among the least reproducibly automated
rules in clinical variant interpretation. The dominant automated implementation counts
ClinVar pathogenic variants within a fixed sequence window, an approach that is
**circular** (ClinVar classifications are themselves ACMG/AMP outputs), window-size
dependent, and undefined on the many genes with sparse ClinVar. We introduce
**AlphaMissense Islands**: contiguous protein segments of consistently high AlphaMissense
pathogenicity, computed proteome-wide and independent of any clinical database. Islands
are obtained from the per-residue mean AlphaMissense score by an optional disorder mask
(AlphaFold pLDDT < 50), a median smoothing filter, and retention of contiguous runs above
the AlphaMissense "ambiguous" boundary (0.564) that are at least 35 residues long. On an
independent benchmark of ClinGen expert-panel curations (ClinGen Evidence Repository,
matched to ClinVar), island-based PM1 reaches **0.70 precision** versus 0.33 (Franklin)
and 0.19 (InterVar) at the best overall F1 (0.59), trading recall for far greater
specificity. A residue-level **error analysis** shows the divergence from expert curation
is systematic and interpretable: two-thirds of false negatives are genuine ClinVar
pathogenic hotspots whose AlphaMissense signal falls just below the island threshold, and
nearly half of false positives are dense pathogenic hotspots that experts simply have not
yet curated. Because islands are orthogonal to ClinVar, they support a **graded PM1**
(Supporting → Moderate → Strong) from cumulative island + in-island ClinVar evidence,
designed to avoid double-counting with PP3. We release the generation script, a
genome-wide BED resource (GRCh38/GRCh37), and **Archipelago**, an interactive browser
covering ~18,400 genes. [CITE:acmg2015][CITE:alphamissense]

---

## 1. Introduction

The ACMG/AMP framework classifies sequence variants by combining evidence criteria of
defined strength. [CITE:acmg2015] **PM1** — a variant "located in a mutational hot spot
and/or critical and well-established functional domain (e.g. active site of an enzyme)
without benign variation" — is defined only qualitatively; the guideline gives no
algorithm for identifying such regions, and automation has been left to individual tool
developers. In practice PM1 is one of the least consistently automated criteria, applied
divergently across classifiers. [CITE:acmg2015][CITE:crosstool2026]

The dominant automated approach infers hotspots from the **density of ClinVar pathogenic
variants** in a sequence window. The clearest documented instance is Franklin's "hotspot
rule": candidate regions are delimited by flanking benign variants and retained only when
they contain pathogenic non-truncating variants, whose density then sets PM1 strength.
[CITE:franklin] This family — implicit in several ACMG engines — has two structural flaws.
First, it is **circular**: ClinVar classifications are themselves products of ACMG/AMP
interpretation, so reusing their density as "independent" PM1 evidence re-uses prior
conclusions. [CITE:clinvar] The circularity is now explicitly recognised: germline
hotspot proposals must guard against pathogenic calls that were themselves informed by
co-located variants. [CITE:haque2025] Second, it is **structurally asymmetric and
fragile** — a fixed window does not track functional boundaries, is sensitive to window
size, collapses where ClinVar is sparse, and gives no benign credit to "coldspots" devoid
of pathogenic variation, causing systematic misclassification (e.g. the BRCA1/BRCA2 coding
coldspots). [CITE:dines2020] A second, more reproducible but coarser family applies PM1
inside **curated functional domains** or hard-coded VCEP hotspot codons (e.g. TP53 codons
175/245/248/249/273/282 plus cancer-hotspot recurrence; PTEN domain rules), which is
gene-specific and does not generalise. [CITE:tp53vcep][CITE:pten][CITE:walsh2018]

Two lines of work aim to break the circularity with signals independent of ClinVar. One
derives sub-genic intolerance from **population data**: missense tolerance ratio (MTR),
regional missense constraint (MPC), constrained coding regions (CCR), homologous missense
constraint at single-residue resolution (HMC), and — most directly — the gnomAD v4
regional-missense-constraint map, which is ClinGen-calibrated so that regions below 20 % of
expected missense variation reach PM1-Moderate support.
[CITE:mtr][CITE:mpc][CITE:ccr][CITE:hmc][CITE:rmc] These signals are powerful but
frequency-based: they require cohorts large enough for depletion to reach significance and
are weakest exactly where clinical genes are rare or short. The other line, well developed
in cancer genomics, calls regional hotspots from the recurrence and 3-D clustering of
observed mutations (cancerhotspots.org, HotMAPS, HotSpot3D).
[CITE:cancerhotspots][CITE:hotmaps][CITE:hotspot3d] Our approach is methodologically closer
to this per-residue-signal → region tradition, but replaces observed recurrence or
population depletion with a *predicted* per-residue effect.

**AlphaMissense** [CITE:alphamissense] provides proteome-wide, per-substitution missense
pathogenicity from a protein language model with structural context — a signal orthogonal
to both clinical databases and allele frequency, and defined even where ClinVar and
population constraint are silent. Its established ACMG use is per-variant **PP3**
(calibrated in-silico evidence). [CITE:pejaver2022] Using AlphaMissense for the *regional*
PM1 criterion is comparatively unexplored, and doing so naïvely risks double-counting with
PP3. We show that a regional summary — an *island* — carries hotspot signal that is
(i) database-independent and non-circular, (ii) transferable to sparse-ClinVar and short
genes, and (iii) separable in principle from the variant's own per-residue score. The
closest prior art is the gnomAD regional-constraint → PM1 preprint [CITE:rmc]: same problem
shape, but a population-depletion signal rather than a predictor, and hence complementary
rather than competing.

**Contributions.** (1) A simple, reproducible definition of AlphaMissense **islands**;
(2) a precomputed genome-wide resource (BED, GRCh38/GRCh37) and its generation script;
(3) **Archipelago**, an open browser over ~18,400 genes; (4) a benchmark against ClinGen
expert curation with a transparent residue-level **error analysis** that makes the
recall/precision trade-off interpretable; and (5) a graded-PM1 framework combining island
membership with in-island ClinVar enrichment, designed to avoid PP3 double-counting.

### 1.1 Related work at a glance

| Approach | Signal | Circular? | Genome-wide? | Representative |
|---|---|---|---|---|
| ClinVar-density window | Submitted P/LP density | **Yes** | Yes | Franklin [CITE:franklin] |
| Curated domain / VCEP codons | Expert annotation | No | No (gene-specific) | TP53, PTEN VCEPs [CITE:tp53vcep][CITE:pten] |
| Population regional constraint | Observed depletion (gnomAD) | No | Yes | RMC→PM1 [CITE:rmc], MTR/MPC/CCR/HMC |
| Cancer hotspot recurrence / 3-D | Observed somatic recurrence | No (somatic) | Yes | cancerhotspots, HotMAPS [CITE:cancerhotspots][CITE:hotmaps] |
| **AlphaMissense islands (this work)** | **Predicted per-residue effect** | **No** | **Yes (~18,400 genes)** | — |

## 2. Methods

### 2.1 Island definition

Per transcript (the AlphaMissense canonical isoform):

1. **Per-residue score** — the mean AlphaMissense pathogenicity across all 19 possible
   substitutions at the residue. [CITE:alphamissense]
2. **Disorder mask (optional)** — residues with AlphaFold pLDDT < 50 set to 0, so islands
   terminate in disordered regions. [CITE:alphafold]
3. **Smoothing** — an edge-aware **median filter, window = 5** residues.
4. **Segmentation** — contiguous runs with smoothed score **> 0.564** (the AlphaMissense
   "ambiguous"/"likely pathogenic" boundary), breaking at the threshold or at gaps in
   residue numbering; runs **< 35 residues** are discarded.

Output is a BED record per island: `chrom, start, end, transcript, uniprot,
island_mean_score, island_length, residue_start, residue_end`. Reference implementation:
`computeAlphaMissenseIslands.py`. Parameter choices (0.564, window 5, min length 35,
pLDDT 50) and a sensitivity sweep are reported in §H1/Supplementary; the released BED was
built with parameters [WITH/WITHOUT pLDDT mask — confirm].

### 2.2 Data & precomputation

Inputs: AlphaMissense (the hg38 per-SNV file for islands and genomic coordinates; the
all-substitutions release for the per-residue saturation heatmap); MANE Select v1.3 for
gene↔transcript mapping [CITE:mane]; the UniProt human reference proteome for domain
features; ClinVar (VEP-annotated, v107, 2025-10) for corroboration [CITE:clinvar]; gnomAD
v2.1.1 regional missense constraint and cancerhotspots.org residue hotspots as optional
context tracks. [CITE:rmc][CITE:cancerhotspots] Per-gene bundles were precomputed for
~18,400 MANE genes and served as static files; each bundle is a compact, gzipped JSON
(quantised heatmap, islands, domains, ClinVar calls, o/e and hotspot tracks).

Because AlphaMissense's canonical transcript differs from MANE for a minority of genes,
ClinVar calls are attached by selecting, per gene, the RefSeq transcript whose reference
amino acids are most concordant with the AlphaMissense sequence, and residues that
disagree are dropped (isoform-concordance gate). The same isoform caveat is surfaced
explicitly in the benchmark (§3.2).

### 2.3 The Archipelago browser

A static web application renders, per gene (Fig. 2): the AlphaMissense saturation heatmap
(all 19 substitutions, continuous ramp), a mean-AM track, islands coloured by ClinVar
overlap (pathogenic-heavy = red, benign-heavy = blue, none = amber), ClinVar P/LP and
B/LB missense lollipops whose **height encodes ClinVar review status** (0–4 gold stars,
from CLNREVSTAT) so expert-panel calls stand out, UniProt domains, gnomAD missense o/e
(regional constraint) and cancer-hotspot tracks, an exon-boundary axis, and summary KPIs
with a ClinVar consequence-breakdown modal. The plot supports drag-to-zoom and PNG export.
A variant page evaluates PM1 for a typed `GENE:p.…`/`GENE:c.…`, and a benchmark explorer
(§3.2) reproduces the PM1 evaluation with a per-variant divergence taxonomy. Live at
archipelago2.vercel.app.

### 2.4 Benchmark & error analysis

**Ground truth.** We use ClinGen expert-panel curations from the ClinGen Evidence
Repository (eRepo), taking PM1 as *applied* when it appears (any strength) in a variant's
"Applied Evidence Codes (Met)". This is an expert, method-independent PM1 label. eRepo
variants are matched to the engine's calls by **ClinVar Variation ID**. [CITE:clingen_svi]

**Prediction.** The engine's PM1 state (`acmg_tags_v2`), which applies PM1 by **genomic
overlap** between the variant and an AlphaMissense island (transcript-agnostic, missense /
in-frame only, no benign downgrade at the binary level). Matching by ClinVar ID means no
protein-isoform mapping is needed for the headline metric.

**Metrics.** True/false positive/negative are computed over the matched set;
precision, recall and F1 are reported for PM1.

**Divergence taxonomy (error analysis).** For every false positive and false negative we
characterise the residue's AlphaMissense neighbourhood from its gene bundle: island
overlap and distance to the nearest island; the count of ClinVar P/LP and B/LB missense
residues within **±10 aa**; the mean AlphaMissense score; and domain membership. Each
false call is assigned one descriptive category (Table 2). Where the eRepo variant maps to
a different isoform than the AlphaMissense canonical sequence (reference amino acid
disagrees), it is labelled *isoform-ambiguous* and left uncategorised rather than guessed.
This analysis is descriptive — it does not re-score any variant.

## 3. Results

### 3.1 PM1 performance vs. existing tools

On the expert-curated benchmark, island-based PM1 is far more **specific** than
window-based classifiers (Table 1): 0.70 precision vs. 0.33 (Franklin) and 0.19
(InterVar), at the best overall F1 (0.59). Window/domain tools reach high recall by
flagging broad regions at a large precision cost. [CITE:intervar]

**Table 1. PM1 against expert curation** *(published SeqOne benchmark; support = 1,270).*

| Tool | Precision | Recall | F1 |
|---|---|---|---|
| **AlphaMissense Islands (SeqOne)** | **0.70** | 0.50 | **0.59** |
| Franklin | 0.33 | 0.86 | 0.47 |
| InterVar | 0.19 | 0.84 | 0.31 |

On the reproducible eRepo/ClinVar-matched set used for the error analysis (10,560 matched
variants; 1,488 with PM1 applied by either side) the engine records 618 TP, 266 FP and
604 FN — precision 0.699, recall 0.506, F1 0.587 — reproducing the headline numbers.
[Add 95 % bootstrap CIs and a paired test vs. comparators — §H4.]

### 3.2 Error analysis: the divergence is systematic and interpretable

Rather than treat every disagreement with experts as noise, we categorise all 266 false
positives and 604 false negatives from their AlphaMissense/ClinVar neighbourhood
(Fig. 3, Table 2).

**False negatives are dominated by sub-threshold hotspots.** Of 604 FN, **403 (67 %)**
are residues surrounded by a genuine ClinVar pathogenic hotspot (≥ 3 P/LP within ±10 aa)
where the AlphaMissense signal was simply not strong or contiguous enough to form an
island — i.e. the recall ceiling is a signal-strength property, not a modelling error.
A further 135 (22 %) have sparse local evidence (little ClinVar, sub-island AM), where
experts likely applied PM1 on non-AM grounds (functional assays, specific-codon VCEP
rules) that a regional predictor cannot see; 29 (5 %) are non-missense (frameshift /
nonsense / splice) to which PM1 — a missense criterion — is by design not applied by the
engine; and 18 (3 %) fall just outside an island boundary.

**False positives are largely defensible hotspots.** Of 266 FP, **121 (45 %)** are
residues in regions dense with ClinVar pathogenic missense (≥ 3 P/LP within ±10 aa): the
island coincides with a real pathogenic hotspot that ClinGen has not (yet) curated PM1
for — a defensible call rather than a clear error. 97 (36 %) are AlphaMissense-only
islands with little ClinVar support (a genuine predictor-vs-panel divergence worth
inspecting), and only 9 (3 %) sit in regions with ≥ 2 nearby benign variants that
weaken the hotspot claim (likely over-calls). 39 FP and 19 FN (58 total) are
isoform-ambiguous and left uncategorised.

**Table 2. Divergence taxonomy of PM1 false calls** *(counts; ±10 aa windows).*

| Class | Category | n | Interpretation |
|---|---|---:|---|
| FP | ClinVar-supported hotspot | 121 | Dense P/LP nearby; ClinGen not-yet-curated → defensible |
| FP | AlphaMissense-only island | 97 | Predictor flags region, little ClinVar support |
| FP | Benign nearby | 9 | ≥2 B/LB nearby; likely over-call |
| FN | Sub-threshold hotspot | 403 | Real ClinVar hotspot, AM below island threshold |
| FN | Sparse evidence | 135 | Non-AM evidence used by the panel |
| FN | Not missense (by design) | 29 | Frameshift/nonsense/splice — PM1 not applicable |
| FN | Just outside an island | 18 | ≤10 aa beyond an island boundary |
| — | Isoform-ambiguous | 58 | Different isoform than AM canonical; uncategorised |

**Recalibrated upper bound.** As an *illustrative* bound (not the benchmark of record),
if the 121 ClinVar-supported-hotspot FPs are counted as correct and the 135 sparse-
evidence + 29 non-missense FNs are treated as out of PM1's scope rather than misses,
precision/recall/F1 move to **0.84 / 0.63 / 0.72**. The point is not the numbers but that
the residual disagreement with experts is concentrated in interpretable, mostly defensible
strata — the profile of a specific, honest criterion rather than a broad, over-firing one.

### 3.3 Orthogonality & illustrative genes

*KIF1A* (Fig. 2) shows motor-domain islands that are pathogenic-enriched (red) while a
C-terminal island reads benign (blue); *PTEN* localises to the phosphatase hotspot; *TP53*
and *BRCA1* (Fig. 4) show tall, star-weighted expert-panel lollipops concentrated in the
DNA-binding and RING/BRCT regions the corresponding VCEPs treat as PM1 territory.
[CITE:tp53vcep][CITE:pten] A formal orthogonality analysis — does island membership add
PM1-relevant signal *beyond the variant's own AlphaMissense/PP3 score* — is required to
defuse double-counting and is reported in §C2 [logistic model: pathogenic ~ variant_AM +
in_island; test in_island after controlling for the point estimate].

### 3.4 Coverage

~18,400 genes are covered. A minority of proteins are absent from AlphaMissense (e.g.
SHANK3/Q9BYB0) or lack a MANE-Select overlap (~1,955 UniProt entries); these are reported
as a limitation.

## 4. Application: grading PM1 by cumulative evidence

Island membership (a structural/predictive hotspot signal) is combined with in-island
ClinVar pathogenic enrichment (empirical corroboration) to modulate PM1 on the ACMG points
scale. [CITE:tavtigian2018][CITE:tavtigian2020]

| Evidence at the residue | PM1 strength |
|---|---|
| Not in an AM island | Not met |
| In an AM island, little/no ClinVar corroboration | Supporting (+1) |
| In an AM island that is pathogenic-enriched (several P/LP, no benign) | Moderate (+2) |
| In a dense, well-established hotspot island (many P/LP, zero benign) | Strong (+4) |
| Island contradicted by benign variation (B/LB-heavy) | Downweight / not met |

The star-weighted ClinVar view (§2.3) makes the corroboration auditable — enrichment
driven by expert-panel (3–4★) submissions is weighted differently from single-submitter
calls. Final cut-offs are gene/VCEP-specific; this grading is presented as a **proposed
framework** pending empirical OddsPath calibration (§C3). [CITE:pejaver2022]

## 5. Discussion

**Advantages.** Islands are non-circular (independent of ClinVar), transferable
(defined on sparse-ClinVar and short genes where population constraint cannot reach
significance), reproducible from a single script + public predictor, and adaptive to
functional boundaries.

**Positioning vs. constraint.** Unlike gnomAD regional missense constraint (frequency-
based), islands are prediction-based and independent of cohort size; the two are
complementary and could be combined. [CITE:rmc] The error analysis shows islands behave as
a high-specificity criterion, with a recall gap that is quantitatively attributable to
sub-threshold AlphaMissense signal.

**Limitations & risks.** PP3/PM1 double-counting (§C2), evidence-strength calibration
(§C3), the benchmark ground-truth definition, parameter robustness (§H1), isoform/coverage
handling, and data licensing.

## 6. Data & code availability

- **Script:** `computeAlphaMissenseIslands.py` (island generation; version-tagged).
- **Resource:** `am_islands_hg38.bed.gz` / `am_islands_hg19.bed.gz` — Zenodo DOI + README
  + column dictionary + exact parameters/data versions.
- **Browser:** Archipelago (open source + precomputed data) — archipelago2.vercel.app;
  archive a tagged release with a DOI.
- **Benchmark:** `benchmark_ACMG/` evaluation harness (reproducible).
- **License:** AlphaMissense predictions are distributed by DeepMind under **CC BY 4.0**
  (the older Zenodo record shows CC BY-NC-SA 4.0); downloading from the CC BY 4.0 source
  permits redistribution of the derived BED and the browser with attribution only.
  [verify exact license string at submission]

## Open items (prioritised)

**Critical** — C1 finalise the ground-truth description (done: ClinGen eRepo expert PM1
calls, matched by ClinVar ID); C2 orthogonality-beyond-PP3 analysis + a co-application
rule; C3 OddsPath calibration or "proposed framework" labelling. **High** — H1 parameter
justification + train/test split + sensitivity sweep (resolve docstring "≥11" vs argparse
"35"); H2 pLDDT-mask ablation + which BED released; H3 comparator method/versions; H4 CIs
+ paired test. **Medium** — hg19 provenance; isoform handling; coverage; RMC agreement;
"islanding vs plain per-residue threshold" ablation; license NOTICE.

## References

*Keys map to `[CITE:key]` above. `[verify]` = confirm ID at submission. Full annotations
and prior-art/novelty verdict in `REFERENCES.md`.*

1. **[acmg2015]** Richards S, et al. ACMG/AMP standards and guidelines. *Genet Med*
   2015;17:405–424. PMID 25741868.
2. **[tavtigian2018]** Tavtigian SV, et al. Bayesian ACMG/AMP framework. *Genet Med*
   2018;20:1054–1060. PMID 29300386.
3. **[tavtigian2020]** Tavtigian SV, et al. Naturally scaled point system. *Hum Mutat*
   2020;41:1734–1737. PMID 32720330.
4. **[pejaver2022]** Pejaver V, et al. ClinGen PP3/BP4 calibration. *Am J Hum Genet*
   2022;109:2163–2177. PMID 36413997.
5. **[alphamissense]** Cheng J, et al. AlphaMissense. *Science* 2023;381:eadg7492.
   PMID 37733863.
6. **[alphafold]** Jumper J, et al. AlphaFold. *Nature* 2021;596:583–589. PMID 34265844.
7. **[alphafold_proteome]** Tunyasuvunakool K, et al. Human-proteome structures. *Nature*
   2021;596:590–596. PMID 34293799.
8. **[mane]** Morales J, et al. MANE. *Nature* 2022;604:310–315. PMID 35388217.
9. **[clinvar]** Landrum MJ, et al. ClinVar. *Nucleic Acids Res* 2020;48:D835–D844.
   PMID 31777943.
10. **[intervar]** Li Q, Wang K. InterVar. *Am J Hum Genet* 2017;100:267–280. PMID 28132688.
11. **[franklin]** Genoox Franklin, "Hotspot rule (PM1)", Help Center (vendor doc).
12. **[varsome]** Kopanos C, et al. VarSome. *Bioinformatics* 2019;35:1978–1980.
    PMID 30376034.
13. **[cardioclassifier]** Whiffin N, et al. CardioClassifier. *Genet Med* 2018;20:1246–1254.
    PMID 29369293.
14. **[charger]** Scott AD, et al. CharGer. *Bioinformatics* 2019;35:865–867. PMID 30102335.
15. **[tapes]** Xavier A, et al. TAPES. *PLoS Comput Biol* 2019;15:e1007453. PMID 31613886.
16. **[pathoman]** Jayakumaran G, et al. PathoMAN. *Genet Med* 2019. PMID 30787465.
17. **[genebe]** Stawiński P, et al. Genebe.net. *Clin Genet* 2024. PMID 38440907.
18. **[genotoscope]** Melidis DP, et al. GenOtoScope. *PLoS Comput Biol* 2022;18:e1009785.
    PMID 36129964.
19. **[crosstool2026]** Ghasemnejad T, et al. Evaluation of ACMG/AMP tools. *Bioinformatics*
    2026;42:btaf623. **[verify PMID]**
20. **[dines2020]** Dines JN, et al. BRCA1/BRCA2 coldspots. *Genet Med* 2020;22:825–830.
    PMID 31911673.
21. **[haque2025]** Haque B, et al. Cancer mutation data for germline classification.
    *PLoS Genet* 2025;21:e1011540. PMID 39761285.
22. **[cancerhotspots]** Chang MT, et al. Recurrent mutations in cancer. *Nat Biotechnol*
    2016;34:155–163. PMID 26619011. *(cancerhotspots.org)*
23. **[cancerhotspots3d]** Chang MT, et al. Functional mutant alleles. *Cancer Discov*
    2018;8:174–183. PMID 29247016.
24. **[hotmaps]** Tokheim C, et al. HotMAPS (3D hotspot regions). *Cancer Res*
    2016;76:3719–3731. PMID 27197156.
25. **[hotspot3d]** Niu B, et al. HotSpot3D. *Nat Genet* 2016;48:827–837. PMID 27294619.
26. **[edriver]** Porta-Pardo E, Godzik A. e-Driver. *Bioinformatics* 2014;30:3109–3114.
    PMID 25064568.
27. **[mutation3d]** Meyer MJ, et al. mutation3D. *Hum Mutat* 2016;37:447–456. PMID 26841357.
28. **[oncokb]** Chakravarty D, et al. OncoKB. *JCO Precis Oncol* 2017;1:1–16. PMID 28890946.
29. **[mtr]** Traynelis J, et al. MTR. *Genome Res* 2017;27:1715–1729. PMID 28864458.
30. **[mpc]** Samocha KE, et al. MPC (regional missense constraint). *bioRxiv* 2017.
    doi:10.1101/148353. *(preprint)*
31. **[ccr]** Havrilla JM, et al. Constrained coding regions. *Nat Genet* 2019;51:88–95.
    PMID 30531870.
32. **[gnomad]** Karczewski KJ, et al. gnomAD constraint. *Nature* 2020;581:434–443.
    PMID 32461654.
33. **[rmc]** Chao KR, et al.; gnomAD. Regional missense intolerance, 730,947 exomes.
    *bioRxiv* 2024. doi:10.1101/2024.04.11.588920. *(preprint; <20% expected missense →
    PM1-Moderate. Closest prior art.)*
34. **[hmc]** Zhang X, et al.; Ware JS. Homologous missense constraint. *Genome Med*
    2024;16:86. PMID 38992748.
35. **[gnocchi]** Chen S, et al. gnomAD v4 genomic constraint (Gnocchi). *Nature*
    2024;625:92–100. PMID 38057664.
36. **[tp53vcep]** Fortuno C, et al. ClinGen TP53 VCEP. *Hum Mutat* 2021;42:223–236.
    PMID 33300245.
37. **[pten]** Mester JL, et al. ClinGen PTEN EP. *Hum Mutat* 2018;39:1581–1592.
    PMID 30311380.
38. **[walsh2018]** Walsh MF, et al. Somatic data for germline classification. *Hum Mutat*
    2018;39:1542–1552. PMID 30311369.
39. **[fayer2021]** Fayer S, et al. Multiplexed functional data (BRCA1/TP53/PTEN).
    *Am J Hum Genet* 2021;108:2248–2258. **[verify PMID]**
40. **[clingen_svi]** ClinGen SVI PM1 guidance. *(No standalone primary doc located —
    VCEP specs + SVI notes; verify at clinicalgenome.org.)*

## Author contributions

Jérôme Audoux, Jon B. [surname/affiliations TBD], et al. — SeqOne.
