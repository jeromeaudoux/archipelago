# Missing pieces to make the AM-Islands paper publishable

Ordered by how badly a reviewer will hit us. ⛔ = likely blocks acceptance until
addressed; 🔶 = expected/strongly advisable; 🟢 = polish.

## ⛔ Critical — the scientific core

### C1. Define the PM1 "ground truth" for the benchmark
The headline (0.70 vs 0.19/0.33) is only meaningful if the reference PM1 label is
defensible and **independent of our method**.
- What set produced support = 1,270? How was "PM1 truly applies at this variant"
  established — ClinGen VCEP expert classifications where PM1 was explicitly invoked
  (from the ClinGen Evidence Repository / eRepo evidence codes)? A held-out ClinVar
  expert-panel set?
- Risk: if the "truth" is itself ClinVar-density-derived, the comparison is unfair in
  our favour. Must state the source, inclusion criteria, and why it's an unbiased PM1
  reference. **This is the first thing a reviewer checks.**

### C2. PP3 vs PM1 double-counting (the killer objection)
ACMG forbids using the same evidence twice. AlphaMissense is already the standard
**per-variant PP3** input. Using AM-derived islands for **PM1** looks like counting the
same signal twice.
- Needed analysis: does island membership carry PM1-relevant signal **independent of
  the variant's own AlphaMissense/PP3 score**? E.g., logistic model: pathogenic ~
  variant_AM_score + in_island; is `in_island` significant after controlling for the
  variant's own score? Or stratify: within variants of equal per-variant AM score, does
  being in an island still enrich for pathogenicity (regional context effect)?
- Decide the framing: (a) islands add orthogonal *regional* evidence beyond the
  per-residue point estimate; and/or (b) a combination rule that prevents
  double-counting (e.g. do not co-apply PP3 and PM1-from-AM, or cap combined strength).
  State this explicitly.

### C3. Evidence-strength calibration (for the graded PM1)
Supporting/Moderate/Strong is currently heuristic. To claim strengths, calibrate to
**local positive predictive value / OddsPath** à la Tavtigian & Pejaver:
- Estimate OddsPath for "variant in island" and for the strata (island + P/LP
  enrichment; dense benign-free island; benign-heavy island) against the truth set, and
  map to Supporting/Moderate/Strong thresholds.
- Until then, phrase the grading as a *proposed framework*, not a validated calibration.

## 🔶 High — expected by reviewers

### H1. Parameter justification + train/test separation
- Justify 0.564 (AM "ambiguous" boundary — cite), window 5, **min length 35**, pLDDT 50.
- ⚠️ Docstring/README inconsistency: module docstring says min length ≥ 11, argparse
  default is 35. Fix and report the value actually used for the released BED.
- Tune parameters on a set **disjoint** from the evaluation set (report the split), with
  a sensitivity sweep (F1 vs threshold/length) so numbers aren't overfit.

### H2. pLDDT-masking ablation + which BED was released
State whether `am_islands_hg38.bed.gz` was built **with or without** pLDDT masking, and
show the effect (island count, PM1 F1) with/without.

### H3. Comparator fairness & reproducibility
- Exactly how InterVar and Franklin assign PM1; tool versions/dates; how Franklin
  (closed-source) calls were obtained; how tool PM1 was matched to the truth label.

### H4. Statistics
95% CIs (bootstrap) on precision/recall/F1; note class balance; the McNemar/paired test
vs comparators on the same variants.

## 🔶 Medium

- **M1. hg19 BED provenance** — native AM hg19 vs liftover (no liftover code in repo);
  state method + any lost regions.
- **M2. Isoform handling** — AM's transcript ≠ MANE for ~28% of genes; islands are keyed
  to the AM transcript. Describe how a variant is mapped to islands at classification
  time (and in the browser, ClinVar NM auto-selected by ref-aa concordance).
- **M3. Coverage** — quantify genes AlphaMissense doesn't cover (e.g. SHANK3 / Q9BYB0)
  and the ~1,955 UniProt entries without a MANE-Select overlap; report as a limitation.
- **M4. Relation to gnomAD regional missense constraint (RMC)** — position islands vs
  RMC (prediction- vs population-based); ideally show agreement/complementarity.
- **M5. Does "islanding" beat a plain per-residue AM threshold?** — ablation: segmentation
  + length filter vs simply "variant's residue mean-AM > 0.564". Justifies the method.
- **M6. License compliance** — AlphaMissense is CC BY-NC-SA 4.0. Derived BED must carry
  share-alike + non-commercial + attribution; confirm the browser (non-commercial) and
  BED redistribution comply; add NOTICE/LICENSE + attribution. (Confirm details from the
  bibliography search.)

## 🟢 Polish / nice-to-have
- G1. 2–3 clinical case studies where island-PM1 changed or clarified a call.
- G2. Pin all data versions (AlphaMissense release, ClinVar v107 2025-10, MANE v1.3,
  gnomAD) + environment for reproducibility.
- G3. Prospective/blind mini-evaluation.
- G4. Figure polish: (i) method schematic (already drafted in the browser Methods page),
  (ii) PM1 benchmark bars, (iii) an example gene panel from Archipelago.

## Community deliverables — to finish before submission
- [ ] **Script** `computeAlphaMissenseIslands.py`: freeze CLI, document params, make
      output deterministic, add a tiny test + the pLDDT-download helper; tag a version.
- [ ] **BED resource**: publish hg38 + hg19 with a README + column dictionary + the exact
      parameters/data versions used; deposit on **Zenodo** for a citable DOI.
- [ ] **Browser (Archipelago)**: source + precomputed data committed; live at
      https://archipelago2.vercel.app; archive a release + DOI.
- [ ] **Benchmark code** (`benchmark_ACMG/`): make the PM1 evaluation reproducible and
      cite it as the validation harness.
- [ ] Author list, affiliations, ORCIDs, funding, competing-interests (SeqOne).

## Suggested minimal path to a submittable short note
1. Nail **C1** (truth set) and **H3** (comparator method) — makes the benchmark credible.
2. Do **C2** (one orthogonality analysis) — defuses the double-counting objection.
3. Add **H1** sensitivity + **H4** CIs.
4. Soften the graded-PM1 to a "proposed framework" unless **C3** calibration is done.
5. Ship script + BED + browser with a Zenodo DOI.
Everything else can be "future work / limitations".
