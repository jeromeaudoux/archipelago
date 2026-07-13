# AlphaMissense Islands: genome-wide regional missense hotspots for non-circular ACMG PM1

*Working draft — Archipelago / SeqOne. Citations marked `[CITE:key]` are resolved in
`REFERENCES.md` (bibliography search in progress). Open items are collected in
`GAPS.md`.*

---

## Title options
1. **AlphaMissense Islands: genome-wide regional missense hotspots for reproducible ACMG PM1**
2. Islands in a sea of missense: turning AlphaMissense into a non-circular PM1 signal
3. A precomputed AlphaMissense hotspot resource and browser for variant classification

## Target format
Short **Application Note** / Brief Communication (≈1,500–2,500 words, 1–2 figures,
1 table). Candidate venues: *Bioinformatics* (Application Note), *Bioinformatics
Advances*, *Genetics in Medicine* (brief report), *GENETICS*/*HGG Advances*; preprint
first on **bioRxiv/medRxiv**. (Venue shortlist pending — see REFERENCES.md §C.)

---

## Abstract (draft)
The ACMG/AMP criterion **PM1** — a variant located in a mutational hotspot or
critical, well-established functional domain — is among the least reproducibly
automated rules. Common implementations count ClinVar pathogenic variants in a fixed
sequence window, which is **circular** (ClinVar classifications are themselves ACMG
outputs) and does not adapt to functional boundaries. We define **AlphaMissense
Islands**: contiguous protein segments of consistently high AlphaMissense
pathogenicity, computed proteome-wide and independent of any clinical database. Islands
are obtained by taking the per-residue mean AlphaMissense score, optionally masking
disordered residues (AlphaFold pLDDT < 50), median-smoothing, and retaining contiguous
runs above the AlphaMissense "ambiguous" boundary (0.564) of ≥ 35 residues. On an
independent benchmark of expert-curated variants, island-based PM1 reaches **0.70
precision** versus 0.19 (InterVar) and 0.33 (Franklin), with the best overall F1 (0.59)
— trading modest recall for far fewer false hotspots. Because islands are orthogonal to
ClinVar, they can be combined with in-island ClinVar enrichment to grade PM1 strength
(Supporting → Moderate → Strong). We release the generation script, a genome-wide BED
resource (GRCh38/GRCh37), and **Archipelago**, an interactive browser covering
~18,400 genes. [CITE:alphamissense][CITE:acmg2015]

---

## 1. Introduction / motivation
- ACMG/AMP 2015 framework; PM1 definition and its known weakness as an automatable
  criterion. [CITE:acmg2015]
- The **circularity problem**: ClinVar-density PM1 uses classifications derived from
  ACMG itself → not independent evidence; also window-size arbitrariness and poor
  transferability to genes with sparse ClinVar. [CITE:clinvar]
- AlphaMissense provides proteome-wide, per-substitution missense pathogenicity from a
  protein language model + structural context — an **orthogonal** signal to clinical
  databases. [CITE:alphamissense]
- Prior use of AlphaMissense in ACMG is mostly **per-variant PP3** (calibrated
  in-silico evidence). [CITE:pejaver2022] Regional/hotspot use is comparatively
  unexplored — *(prior-art check pending, REFERENCES.md §B)*. Related prior concept:
  gnomAD **regional missense constraint**, which is population- rather than
  prediction-based. [CITE:rmc][CITE:gnomad]
- Contribution: (i) a simple, reproducible definition of AlphaMissense **islands**;
  (ii) a genome-wide precomputed resource + generation script; (iii) an open browser;
  (iv) a benchmark showing large precision gains for PM1; (v) a scheme to grade PM1
  strength by cumulative island + ClinVar evidence.

## 2. Methods
### 2.1 Island definition (the algorithm)
Per transcript (AlphaMissense's canonical isoform):
1. **Per-residue score** = mean AlphaMissense pathogenicity over all substitutions at
   the residue. [CITE:alphamissense]
2. **Disorder mask (optional)**: residues with AlphaFold pLDDT < 50 set to 0, so
   islands terminate in disordered regions. [CITE:alphafold]
3. **Smoothing**: edge-aware **median filter, window = 5** residues.
4. **Segmentation**: contiguous runs with smoothed score **> 0.564** (AlphaMissense
   "ambiguous" boundary), breaking at the threshold or at gaps in residue numbering;
   runs **< 35 residues** discarded.

Output = BED (chrom, start, end, transcript, uniprot, island_mean_score,
island_length, residue_start, residue_end). Reference implementation:
`computeAlphaMissenseIslands.py`. *(Parameter justification / sensitivity — GAPS.md.)*

### 2.2 Data & precomputation
- Inputs: AlphaMissense hg38 (per-SNV, for islands + coordinates) and the all-19
  aa_substitutions release (for the per-residue saturation heatmap); MANE Select for
  symbol↔transcript mapping; UniProt human proteome for domains; ClinVar (VEP-annotated,
  v107 2025-10) for corroboration. [CITE:alphamissense][CITE:mane][CITE:clinvar]
- Per-gene bundles precomputed for ~18,400 MANE genes; served statically.

### 2.3 The Archipelago browser
Static web app: per-gene AlphaMissense saturation heatmap, mean-AM track, islands
(colored by ClinVar overlap), ClinVar P/LP & B/LB missense lollipops, UniProt domains,
and KPIs. URL: https://archipelago2.vercel.app.

### 2.4 Benchmark
- Truth set: expert-curated variants *(exact source + how PM1 ground-truth is defined —
  GAPS.md, critical)*; comparators InterVar and Franklin. [CITE:intervar]
- Metric: per-criterion precision/recall/F1 for PM1 against the reference PM1 assignment.

## 3. Results
### 3.1 PM1 performance
| Tool | Precision | Recall | F1 |
|---|---|---|---|
| **AlphaMissense Islands (SeqOne)** | **0.70** | 0.50 | **0.59** |
| Franklin | 0.33 | 0.86 | 0.47 |
| InterVar | 0.19 | 0.84 | 0.31 |

(support = 1,270). Window-based tools achieve high recall by flagging broad regions,
at a large precision cost; islands are far more specific. *(Add 95% CIs — GAPS.md.)*

### 3.2 Orthogonality & illustrative genes
- KIF1A (motor-domain islands, red/pathogenic; C-terminal benign), PTEN (phosphatase
  hotspot), TP53 — screenshots/figure from the browser.
- Argument that islands are not redundant with per-variant PP3 *(needs quantification —
  GAPS.md, the double-counting concern).*

### 3.3 Coverage
~18,400 genes; note genes AlphaMissense does not cover (e.g. some isoforms/proteins).

## 4. Application: grading PM1 by cumulative evidence
Combine island membership (structural/predictive hotspot) with **in-island ClinVar
pathogenic enrichment** (empirical corroboration) to modulate PM1 on the ACMG points
scale (Supporting = 1, Moderate = 2, Strong = 4). [CITE:tavtigian2018]
- Island only → Supporting; island + P/LP enrichment (no benign) → Moderate; dense,
  benign-free hotspot → Strong; benign-heavy island → downweight.
- Final cut-offs are gene/VCEP-specific. [CITE:vcep] *(Needs empirical calibration —
  GAPS.md.)*

## 5. Discussion
- Advantages: non-circular, transferable (incl. sparse-ClinVar genes), reproducible,
  adapts to functional boundaries.
- Limitations & risks (see GAPS.md): PP3/PM1 double-counting, strength calibration,
  benchmark ground-truth definition, parameter robustness, isoform/coverage, license.

## 6. Data & code availability
- **Script**: `computeAlphaMissenseIslands.py` (island generation).
- **Resource**: `am_islands_hg38.bed.gz` / `am_islands_hg19.bed.gz`.
- **Browser**: Archipelago (source + precomputed data), https://archipelago2.vercel.app.
- License note: AlphaMissense is **CC BY-NC-SA 4.0** — implications for the derived BED
  and browser under review *(REFERENCES.md, license note)*.

## Author contributions
Jérôme Audoux, Jon B. [surname/affiliations TBD], et al. — SeqOne.
