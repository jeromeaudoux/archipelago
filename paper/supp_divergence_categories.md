# Supplementary Table S2 — Divergence-category definitions (PM1 error analysis §3.3)

Every PM1 false positive (FP) and false negative (FN) on the eRepo/ClinVar-matched set is
assigned exactly one descriptive category from its AlphaMissense neighbourhood: island
overlap and distance to the nearest island, ClinVar P/LP and B/LB missense-residue counts
within ±10 aa, mean AlphaMissense score, and domain membership. The assignment is
descriptive and does not re-score any variant. Categories and their live definitions are
browsable in the PM1 benchmark explorer (archipelago2.vercel.app, "Where do the engine and
ClinGen diverge?").

| Class | Category | n | Definition (residue neighbourhood) |
|---|---|---:|---|
| FP | ClinVar-supported hotspot | 121 | In an island; ≥3 ClinVar P/LP missense within ±10 aa and <2 B/LB — a real pathogenic hotspot ClinGen has not (yet) curated PM1 for (defensible). |
| FP | AlphaMissense-only island | 97 | In an island; <3 P/LP nearby — predictor flags the region with little ClinVar support (genuine predictor-vs-panel divergence). |
| FP | Benign nearby | 9 | In an island; ≥2 B/LB missense within ±10 aa — benign variation weakens the hotspot claim (likely over-call). |
| FN | Sub-threshold hotspot | 403 | Not in an island but ≥3 P/LP within ±10 aa — a real ClinVar pathogenic hotspot whose AlphaMissense signal did not clear the island threshold. |
| FN | Just outside an island | 18 | Not in an island but ≤10 aa beyond an island boundary (boundary sensitivity). |
| FN | Sparse evidence | 135 | Not in an island, <3 P/LP nearby, sub-island AM — the panel likely used non-AM evidence (functional assays, specific-codon VCEP rules). |
| FN | Not missense (by design) | 29 | Frameshift / nonsense / splice / non-single-residue — PM1 is a missense criterion, so PM1_AM is correctly not applied. |
| — | Isoform-ambiguous | 58 | eRepo variant maps to a different isoform than the AlphaMissense canonical sequence (reference amino acid disagrees); left uncategorised rather than guessed (39 FP, 19 FN). |

*Counts are for the ClinVar-matched eRepo set (266 FP, 604 FN). Windows: ±10 aa;
"hotspot" threshold ≥3 P/LP. Generated alongside `public/data/pm1_benchmark.json`.*
