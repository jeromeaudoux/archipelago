# Do islands capture *critical functional regions* (PM1) or generic conservation?

*Reproduce: `python3 explore/island_functional_enrichment.py`. Figure:
`paper/figures/fig_functional_enrichment.svg` / `.png`.*

## Question

PM1 is *"a mutational hot spot and/or critical and well-established functional domain."* Do
AlphaMissense islands actually land on functional regions, or are they just any conserved
stretch a predictor lights up? We test island vs non-island residues, **within the same
island-bearing proteins**, for three markers of functional criticality that are each computed
**independently of AlphaMissense**.

## Result (6,148 island-bearing proteins; island vs non-island residues)

| Independent annotation | In island | In non-island | Enrichment OR [95% CI] |
|---|--:|--:|--:|
| UniProt functional domain | 33 % | 17 % | **2.3× [2.1–2.5]** |
| Cancer hotspot (±2 aa, cancerhotspots.org) | 0.03 % | 0.01 % | **4.3× [2.4–7.6]** |
| gnomAD regional missense constraint (o/e < 0.6) | 43 % | 14 % | **4.7× [4.4–4.9]** |

CIs are gene-cluster bootstraps.

## Interpretation

Island residues are **2–5× enriched** for curated functional domains, recurrent somatic
hotspots, and population-depleted regions — three orthogonal, non-AlphaMissense annotations
of functional importance. Convergent enrichment across all three argues that islands mark
genuine **critical functional regions**, i.e. exactly what PM1 names, rather than the
predictor's own conservation signal alone. (The cancer-hotspot prevalence is tiny in absolute
terms — hotspots are rare — but the enrichment is clear.)

## Why this supports PM1 *specifically* (and not other ACMG criteria)

The evidence is **region-level**, matching PM1's granularity. It does not map onto:
- **PP2** — a *gene-level* property (missense is a common disease mechanism + low benign
  missense rate); islands are sub-genic and say nothing about a gene's mechanism.
- **PP3/BP4** — a *per-variant* predictor; islands are substitution-agnostic regions (and §3.4
  shows they are non-redundant with the deployed PP3).
- **PM5/PS1** — *residue/allele*-specific; islands are contiguous regions, not single codons.

So the mapping to PM1 is structural: the unit of the signal (a critical functional region) is
the unit of the criterion. Islands operationalise the *"critical, well-established functional
domain"* limb of PM1; the narrow recurrent-hotspot limb (which islands under-capture) is
exactly the source of the sub-threshold-hotspot false negatives in §3.3.

## Caveats

Within-protein comparison controls for protein-level composition but not for the fact that
AlphaMissense and gnomAD constraint both reflect intolerance (they are different data sources
— prediction vs population — but correlated); UniProt domains and cancer hotspots are fully
independent. "UniProt functional domain" is a broad annotation; finer active-/binding-site
enrichment would sharpen the "critical" claim further if a site-level UniProt feature file is
added.

## Paper use

Results §3.6 (Fig. 7) + the "Specific to PM1" Discussion argument.
