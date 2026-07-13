# References, prior art & venues

*Verified via primary sources; items marked **[verify]** still need DOI/PMID confirmed at
submission.*

## A. Bibliography
1. **AlphaMissense** — Cheng J, et al. *Accurate proteome-wide missense variant effect
   prediction with AlphaMissense.* Science 2023;381(6664):eadg7492.
   doi:10.1126/science.adg7492 · PMID 37733863. *(Do NOT cite adj8672 — that's the
   Perspective.)*
2. **AlphaFold / pLDDT** — Jumper J, et al. Nature 2021;596:583–589.
   doi:10.1038/s41586-021-03819-2 · PMID 34265844.
3. **MANE Select** — Morales J, et al. Nature 2022;604:310–315.
   doi:10.1038/s41586-022-04558-8 · PMID 35388217.
4. **ACMG/AMP 2015** — Richards S, et al. Genet Med 2015;17(5):405–424.
   doi:10.1038/gim.2015.30 · PMID 25741868.
5. **Bayesian ACMG** — Tavtigian SV, et al. Genet Med 2018;20(9):1054–1060.
   doi:10.1038/gim.2017.210 · PMID 29300386.
6. **ACMG points** — Tavtigian SV, et al. Hum Mutat 2020;41(10):1734–1737.
   doi:10.1002/humu.24088. **[verify PMID]**
7. **PP3/BP4 calibration (incl. AlphaMissense)** — Pejaver V, et al. Am J Hum Genet
   2022;109(12):2163–2177. doi:10.1016/j.ajhg.2022.10.013 · PMID 36413997.
   *(Follow-up: Schmidt A, et al. GIM 2025, doi 10.1016/j.gim.2025.101418 [verify];
   bioRxiv 2024.09.17.611902.)*
8. **ClinVar** — Landrum MJ, et al. Nucleic Acids Res 2020;48(D1):D835–D844.
   doi:10.1093/nar/gkz972 · PMID 31777943. *(also 2018: doi 10.1093/nar/gkx1153.)*
9. **gnomAD constraint (v2)** — Karczewski KJ, et al. Nature 2020;581:434–443.
   doi:10.1038/s41586-020-2308-7 · PMID 32461654.
10. **Regional missense constraint** — Samocha KE, et al. bioRxiv 2017,
    doi:10.1101/148353 (MPC). gnomAD v4 RMC (730,947 exomes) bioRxiv 2024,
    doi:10.1101/2024.04.11.588920 — *maps regional o/e depletion to calibrated PM1
    strengths; the key work to contrast.*
11. **gnomAD v4 genomic constraint (Gnocchi)** — Chen S, et al. Nature 2024;625:92–100.
    doi:10.1038/s41586-023-06045-0. **[verify PMID ~38057664]**
12. **InterVar** — Li Q, Wang K. Am J Hum Genet 2017;100(2):267–280.
    doi:10.1016/j.ajhg.2017.01.004 · PMID 28132688.
13. **VCEP PM1 examples** — TP53: Fortuno C, et al. Hum Mutat 2021;42(3):223–236
    (doi:10.1002/humu.24152). PTEN: Mester JL, et al. Hum Mutat 2018;39(11):1581–1592
    (doi:10.1002/humu.23636).
14. **Cross-tool ACMG/PM1 benchmarks (motivation)** — Bioinformatics 2025,
    doi:10.1093/bioinformatics/btaf623; Brief Bioinform 2025, bbaf545.

## B. Prior art & novelty
- **Biggest overlap:** gnomAD RMC→PM1 preprint (#10). Concept overlaps (sub-genic
  regions → graded PM1) but the signal is **population depletion (o/e)**, not an
  in-silico predictor. **Differentiate explicitly:** AlphaMissense islands are
  predictor-derived (work where cohorts are too small for constraint to reach
  significance) and orthogonal to allele frequency.
- AlphaMissense in ACMG elsewhere is **per-variant PP3 only** (e.g. Kurtovic-Kozaric
  et al., Front Genet 2024, doi:10.3389/fgene.2024.1487608). No genome-wide AM *region*
  resource for PM1 exists.
- Adjacent: AlphaMissenseR (Bioinform Adv 2025, vbaf093) — access/viz, no PM1; MTR3D,
  structure-based intolerance — not PM1, not a BED/PM1 resource.
- **Verdict:** the combination — precomputed genome-wide AM islands + graded PM1 + open
  browser + BED, designed non-circular with PP3 — is **unpublished**. Foreground the
  AM-vs-constraint distinction and the non-circularity design.

## C. License (important, changed)
- DeepMind's official AlphaMissense repo now states predictions are **CC BY 4.0** (code
  Apache-2.0); Ensembl VEP plugin header moved from CC BY-NC-SA 4.0 (rel 110) to CC BY
  4.0 (rel 115). The **old Zenodo record still shows CC BY-NC-SA 4.0**.
- **Recommendation:** download from the DeepMind CC BY 4.0 source and cite host + license
  + access date. Under CC BY 4.0 the derived BED + browser need **attribution only** (no
  NC, no SA). **[verify the exact license string on your download at submission.]**

## D. Venues
Priority preprint on **bioRxiv/medRxiv** (protect against the RMC-PM1 preprint), then:
- Tooling lead → **Bioinformatics Advances** (OA) or **Bioinformatics** Application Note.
- Resource/browser lead → **NAR Web Server/Database** (max reach) or **GigaScience/Database**.
- Clinical lead → **HGG Advances** or **GIM Open**.
