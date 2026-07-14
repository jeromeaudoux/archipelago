# References, prior art & venues

*Verified against primary sources (PubMed / publisher / PMC) during the prior-work search.
Items marked **[verify]** still need a DOI/PMID confirmed at submission. Citation keys
match `[CITE:key]` in `DRAFT-v2.md`.*

## A. Bibliography

### Framework & evidence strength
- **[acmg2015]** Richards S, Aziz N, Bale S, et al. Standards and guidelines for the
  interpretation of sequence variants (ACMG/AMP). *Genet Med* 2015;17(5):405–424.
  doi:10.1038/gim.2015.30 · PMID 25741868. *(PM1 defined qualitatively; no algorithm — the
  root of the automation gap.)*
- **[tavtigian2018]** Tavtigian SV, et al. Modeling the ACMG/AMP guidelines as a Bayesian
  classification framework. *Genet Med* 2018;20(9):1054–1060. doi:10.1038/gim.2017.210 ·
  PMID 29300386.
- **[tavtigian2020]** Tavtigian SV, et al. Fitting a naturally scaled point system to the
  ACMG/AMP guidelines. *Hum Mutat* 2020;41(10):1734–1737. doi:10.1002/humu.24088 ·
  PMID 32720330.
- **[pejaver2022]** Pejaver V, et al.; ClinGen SVI. Calibration of computational tools for
  missense variant pathogenicity and ClinGen recommendations for PP3/BP4. *Am J Hum Genet*
  2022;109(12):2163–2177. doi:10.1016/j.ajhg.2022.10.013 · PMID 36413997. *(Local-PPV
  calibration = the template for grading PM1 strength; also the PP3 use to distinguish
  from.)* Follow-up (more tools): AJHG 2024, PMID 39345488 *(check if AlphaMissense
  recalibrated for PP3/BP4).*
- **[pvs1]** Abou Tayoun AN, et al.; ClinGen SVI. Recommendations for interpreting the
  loss-of-function PVS1 ACMG/AMP criterion. *Hum Mutat* 2018;39(11):1517–1524.
  doi:10.1002/humu.23626 · PMID 30192042. *(Example of SVI operationalising a hard
  criterion — the contrast being that no such PM1 specification exists.)*
- **[clingen_svi]** ClinGen Sequence Variant Interpretation Working Group PM1 guidance —
  *no standalone primary document located; VCEP specs + secondary SVI notes only.*
  **[verify at clinicalgenome.org — e.g. PM1+PP3 combination cap; PM1↓ when PS1/PM5 apply.]**

### Predictors
- **[alphamissense]** Cheng J, Novati G, Pan J, et al. Accurate proteome-wide missense
  variant effect prediction with AlphaMissense. *Science* 2023;381(6664):eadg7492.
  doi:10.1126/science.adg7492 · PMID 37733863. *(Research article — NOT the Perspective
  adj8672.)*
- **[alphafold]** Jumper J, et al. Highly accurate protein structure prediction with
  AlphaFold. *Nature* 2021;596:583–589. doi:10.1038/s41586-021-03819-2 · PMID 34265844.
- **[alphafold_proteome]** Tunyasuvunakool K, et al. Highly accurate protein structure
  prediction for the human proteome. *Nature* 2021;596:590–596.
  doi:10.1038/s41586-021-03828-1 · PMID 34293799. *(The proteome-scale companion; cite
  whichever pLDDT source was actually used for the disorder mask.)*
- **[mane]** Morales J, et al. A joint NCBI/EMBL-EBI transcript set (MANE). *Nature*
  2022;604:310–315. doi:10.1038/s41586-022-04558-8 · PMID 35388217.
- **[clinvar]** Landrum MJ, et al. ClinVar: improvements to accessing data. *Nucleic Acids
  Res* 2020;48(D1):D835–D844. doi:10.1093/nar/gkz972 · PMID 31777943.

### Automated ACMG classifiers (and their PM1 rule)
- **[intervar]** Li Q, Wang K. InterVar: clinical interpretation by the 2015 ACMG/AMP
  guidelines. *Am J Hum Genet* 2017;100(2):267–280. doi:10.1016/j.ajhg.2017.01.004 ·
  PMID 28132688. *(Automates 18/28 criteria; exact PM1 region source needs code-level
  confirmation — reported to over-flag PM1 vs eRepo.)*
- **[franklin]** Genoox Franklin — "Hotspot rule (PM1)", Help Center (vendor documentation;
  no peer-reviewed primary paper). help.genoox.com/en/articles/4230441. *(Region delimited
  by flanking benign variants, retained if it contains pathogenic non-truncating variants,
  density → PM1 strength — the canonical circular ClinVar-density method to contrast.)*
- **[varsome]** Kopanos C, et al. VarSome. *Bioinformatics* 2019;35(11):1978–1980.
  doi:10.1093/bioinformatics/bty897 · PMID 30376034. *(PM1 rule not detailed in the paper.)*
- **[cardioclassifier]** Whiffin N, et al. CardioClassifier. *Genet Med*
  2018;20(10):1246–1254. doi:10.1038/gim.2017.258 · PMID 29369293. *(Disease-specific rule
  engineering — contrast to genome-wide generic PM1.)*
- **[charger]** Scott AD, et al. CharGer. *Bioinformatics* 2019;35(5):865–867.
  PMID 30102335.
- **[tapes]** Xavier A, et al. TAPES. *PLoS Comput Biol* 2019;15(10):e1007453.
  PMID 31613886.
- **[pathoman]** Jayakumaran G, et al. PathoMAN. *Genet Med* 2019.
  doi:10.1038/s41436-019-0463-8 · PMID 30787465.
- **[genebe]** Stawiński P, et al. Genebe.net. *Clin Genet* 2024. doi:10.1111/cge.14516 ·
  PMID 38440907.
- **[genotoscope]** Melidis DP, et al. GenOtoScope (hearing loss). *PLoS Comput Biol*
  2022;18(9):e1009785. PMID 36129964. *(Implements PM1 via curated hotspot lists.)*
- **[crosstool2026]** Ghasemnejad T, et al. Comprehensive evaluation of ACMG/AMP-based
  variant classification tools. *Bioinformatics* 2026;42(2):btaf623.
  doi:10.1093/bioinformatics/btaf623. **[verify PMID]** *(Benchmarks Franklin, InterVar,
  TAPES, Genebe. The "InterVar over-flags PM1" claim came via a summary layer — verify
  against the paper.)* *(A separate 2025 Brief Bioinform PM1 paper was NOT confirmed to
  exist — drop unless found.)*

### Circularity & hotspot-asymmetry critiques
- **[dines2020]** Dines JN, et al. Systematic misclassification of missense variants in
  BRCA1/BRCA2 "coldspots". *Genet Med* 2020;22(5):825–830.
  doi:10.1038/s41436-019-0740-6 · PMID 31911673. *(Hotspots credited, coldspots not —
  motivates a symmetric, genome-wide signal.)*
- **[haque2025]** Haque B, et al. Leveraging cancer mutation data to inform pathogenicity
  classification of germline missense variants. *PLoS Genet* 2025;21(1):e1011540.
  doi:10.1371/journal.pgen.1011540 · PMID 39761285. *(Explicitly names ClinVar circularity;
  extends Walsh 2018 to Mendelian genes.)*

### Mutational hotspots (cancer; transferable concept)
- **[cancerhotspots]** Chang MT, et al. Identifying recurrent mutations in cancer.
  *Nat Biotechnol* 2016;34(2):155–163. doi:10.1038/nbt.3391 · PMID 26619011.
- **[cancerhotspots3d]** Chang MT, et al. Accelerating discovery of functional mutant
  alleles. *Cancer Discov* 2018;8(2):174–183. doi:10.1158/2159-8290.CD-17-0321 ·
  PMID 29247016.
- **[hotmaps]** Tokheim C, et al. Exome-scale discovery of hotspot regions using 3D protein
  structure (HotMAPS). *Cancer Res* 2016;76(13):3719–3731.
  doi:10.1158/0008-5472.CAN-15-3190 · PMID 27197156.
- **[hotspot3d]** Niu B, et al. Protein-structure-guided discovery of functional mutations
  (HotSpot3D). *Nat Genet* 2016;48(8):827–837. doi:10.1038/ng.3586 · PMID 27294619.
- **[edriver]** Porta-Pardo E, Godzik A. e-Driver. *Bioinformatics* 2014;30(21):3109–3114.
  doi:10.1093/bioinformatics/btu499 · PMID 25064568.
- **[mutation3d]** Meyer MJ, et al. mutation3D. *Hum Mutat* 2016;37(5):447–456.
  doi:10.1002/humu.22963 · PMID 26841357.
- **[oncokb]** Chakravarty D, et al. OncoKB. *JCO Precis Oncol* 2017;1:1–16.
  doi:10.1200/PO.17.00011 · PMID 28890946.

### Regional / sub-genic constraint (closest neighbours)
- **[mtr]** Traynelis J, et al. Gene-customized missense interpretation (MTR).
  *Genome Res* 2017;27(10):1715–1729. doi:10.1101/gr.226589.117 · PMID 28864458.
- **[mpc]** Samocha KE, et al. Regional missense constraint improves deleteriousness
  prediction (MPC). *bioRxiv* 2017. doi:10.1101/148353. *(Preprint only — confirm no journal
  version.)*
- **[ccr]** Havrilla JM, et al. A map of constrained coding regions (CCR). *Nat Genet*
  2019;51(1):88–95. doi:10.1038/s41588-018-0294-6 · PMID 30531870.
- **[gnomad]** Karczewski KJ, et al. Mutational constraint spectrum (gnomAD, 141,456
  humans). *Nature* 2020;581:434–443. doi:10.1038/s41586-020-2308-7 · PMID 32461654.
- **[rmc]** Chao KR, Wang L, Panchal R, et al.; gnomAD. The landscape of regional missense
  mutational intolerance quantified from 730,947 exomes. *bioRxiv* 2024.
  doi:10.1101/2024.04.11.588920 (PMC11030311). **Preprint — not yet peer-reviewed.**
  *(THE work to contrast: ClinGen-calibrated, <20% expected missense → PM1-Moderate. Signal
  = population depletion vs. our predictor; complementary.)*
- **[hmc]** Zhang X, Theotokis PI, ..., Ware JS. Genetic constraint at single-amino-acid
  resolution in protein domains (Homologous Missense Constraint). *Genome Med* 2024;16:86.
  doi:10.1186/s13073-024-... · PMID 38992748 (medRxiv 2022 doi:10.1101/2022.02.16.22271023).
  *(Novelty-adjacent: single-residue-resolution constraint for missense prioritisation —
  read full text before finalising novelty claims.)*
- **[gnocchi]** Chen S, et al. A genomic mutational constraint map (gnomAD v4, Gnocchi).
  *Nature* 2024;625:92–100. doi:10.1038/s41586-023-06045-0 · PMID 38057664. *(Non-coding /
  genomic constraint — cite for completeness, do not conflate with RMC.)*

### VCEP / domain-based PM1 & functional data
- **[tp53vcep]** Fortuno C, et al.; ClinGen TP53 VCEP. *Hum Mutat* 2021;42(3):223–236.
  doi:10.1002/humu.24152 · PMID 33300245. *(PM1 = specific hotspot codons + cancerhotspots
  ≥10 occurrences; VUS 28%→12%.)*
- **[pten]** Mester JL, et al.; ClinGen PTEN EP. *Hum Mutat* 2018;39(11):1581–1592.
  doi:10.1002/humu.23636 · PMID 30311380.
- **[walsh2018]** Walsh MF, et al. Integrating somatic variant data for germline
  classification in cancer-predisposition genes. *Hum Mutat* 2018;39(11):1542–1552.
  doi:10.1002/humu.23640 · PMID 30311369. *(Origin of cancerhotspots-informed germline PM1.)*
- **[rasopathy2025]** ClinGen RASopathy VCEP. Updated ACMG/AMP specifications for variant
  interpretation and gene curations from the ClinGen RASopathy expert panels. *Genet Med
  Open* 2025 (PMC12151217). **[verify DOI/PMID]** *(Evidence the field is disambiguating PM1:
  recommends reserving PM1 for well-established functional domains and using PM5_Strong for
  mutational hot-spots — directly relevant to how the revised guidelines may frame PM1.)*
- **[fayer2021]** Fayer S, et al. Systematic integration of multiplexed functional data
  (BRCA1/TP53/PTEN). *Am J Hum Genet* 2021;108(12):2248–2258.
  doi:10.1016/j.ajhg.2021.11.001 · PMID **[verify]**. *(MAVE route to regional intolerance;
  contrast: functional ground truth, gene-by-gene, not genome-wide.)*

## B. Prior art & novelty verdict

- **Closest prior art:** the gnomAD RMC→PM1 preprint **[rmc]** (2024). Same problem shape
  (sub-genic region → calibrated PM1), different signal (population depletion vs.
  AlphaMissense prediction). **Differentiate explicitly and frame as complementary** —
  islands are defined where cohorts are too small for constraint significance.
- **Novelty-adjacent:** HMC **[hmc]** (single-residue constraint, published Genome Med
  2024) — read full text; different substrate (homology-aggregated population constraint).
- **Methodological lineage:** cancer 3-D hotspot clustering **[hotmaps][hotspot3d]** shows
  "per-residue signal → region" is well-trodden, but not for germline PM1 from a predictor.
- **Direct comparator (product):** Franklin's documented ClinVar-density hotspot rule
  **[franklin]** — the circular method we position against (vendor doc, not literature).
- **Verdict:** the specific combination — precomputed genome-wide **AlphaMissense** islands
  + graded PM1 + open browser + BED, designed non-circular and separable from PP3 — is
  unpublished. Foreground (i) predictor-vs-population/ClinVar independence, and (ii) the
  honest residue-level error analysis.

## C. Open citation checks
- `[verify]`: crosstool2026 PMID; fayer2021 PMID; hmc exact Genome Med DOI; MPC/RMC journal
  status; a real ClinGen-SVI PM1 primary doc; whether a 2025 Brief Bioinform PM1 paper
  exists.

## D. Venues
Preprint on **bioRxiv/medRxiv** first (priority vs. [rmc]), then:
- Tooling → **Bioinformatics Advances** (OA) / **Bioinformatics** Application Note.
- Resource/browser → **NAR Web Server/Database** / **GigaScience/Database**.
- Clinical → **HGG Advances** / **GIM Open**.

## E. License
AlphaMissense predictions: DeepMind repo/GCS state **CC BY 4.0** (code Apache-2.0); older
Zenodo record shows CC BY-NC-SA 4.0. Download from the CC BY 4.0 source → derived BED +
browser need **attribution only**. **[verify the exact license string on the download at
submission.]**
