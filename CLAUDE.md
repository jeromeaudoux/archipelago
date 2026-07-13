# am-islands-viewer — project memory

## What this is

A **standalone, Vercel-deployable web app**: a user types a human gene symbol and gets an
interactive viewer showing, aligned on the protein sequence:

- an **AlphaMissense missense saturation heatmap** (20 amino acids × every residue,
  colored by continuous AM pathogenicity 0–1),
- a **mean-AM** per-residue track,
- **AlphaMissense islands** (contiguous high-pathogenicity regions), colored by their
  **ClinVar overlap** (red = pathogenic-heavy, blue = benign-heavy, amber = no ClinVar),
- **ClinVar** Pathogenic/Likely-pathogenic (lollipops up) and Benign/Likely-benign (down),
- **UniProt protein domains**.

Origin: generalized from a one-off KIF1A artifact. The goal is coverage of **all ~19k
human MANE-Select genes**.

## Core design constraint

Raw sources are huge (AlphaMissense hg38 ≈ 624 MB gz; ClinVar VEP VCF ≈ 679 MB gz). They
**never ship to the browser**. Instead we **precompute one compact gzipped JSON bundle per
gene** offline and commit them to `public/data/genes/<SYMBOL>.json.gz` (~few KB each,
~60–90 MB total). The browser downloads only the bundle for the gene being viewed and
inflates it client-side (`DecompressionStream`).

## Raw data sources (NOT committed — live outside the repo)

| Source | Path | Role |
|---|---|---|
| AlphaMissense **aa_substitutions** | `~/Downloads/AlphaMissense_aa_substitutions.tsv.gz` | **heatmap source** — all 19 substitutions per residue, keyed by UniProt only (no coords/transcript) |
| AlphaMissense hg38 | `~/Downloads/AlphaMissense_hg38.tsv.gz` | used only to derive `UniProt→symbol/ensembl` (genomic overlap w/ MANE); the SNV-based source of the islands bed |
| AM islands | `../acmgscore-v2/data/am_islands_hg38.bed.gz` | precomputed islands, joined by UniProt (`uniprot_id` col) |
| Protein domains | `../acmgscore-v2/UP000005640_9606_domain.bed` | UniProt features (`X<start>-Y<end>; NAME` in last col), keyed by UniProt |
| ClinVar VEP | `../Data/resources/human/clinvar-vep/GRCh38/107_2025-10/clinvar.vep.vcf.gz` | CSQ: SYMBOL, Feature (RefSeq NM), Protein_position, Amino_acids; CLNSIG on record |
| MANE summary | `../Data/MANE.GRCh38.v1.3.summary.txt.gz` | mapping backbone: symbol ↔ Ensembl(ENST) ↔ RefSeq(NM) ↔ UniProt |

Paths are configurable at the top of `scripts/config.py`.

AM file columns: `CHROM POS REF ALT genome uniprot_id transcript_id protein_variant
am_pathogenicity am_class`. `protein_variant` e.g. `V2L` (ref+pos+alt).

## Why MANE

AlphaMissense is keyed by Ensembl transcript (one canonical per protein) and carries no
gene symbol; ClinVar VEP is keyed by RefSeq NM. The MANE Select summary bridges
`Ensembl ↔ symbol ↔ RefSeq`, so we pick exactly **one transcript per gene** and can join
all three sources. **Always join on unversioned ENST/NM bases** (versions drift between
releases); keep versioned IDs only as display metadata.

## Bundle schema (`public/data/genes/<SYMBOL>.json.gz`)

```jsonc
{
  "gene": "KIF1A", "uniprot": "Q12756",
  "ensembl": "ENST00000650053.1", "refseq": "NM_004321.8",
  "length": 1690,
  "aa_order": "HKRDESTNQCGPAVILMFYW",      // fixed heatmap row order
  "mean": [0.72, ...],                      // per-position mean AM (len N; null if unscored)
  "ref":  "MAGA...",                         // reference aa string (len N)
  "grid": "<base64 of N*20 bytes>",          // per-substitution score, row-major by aa_order
  // byte: 0..250 => score/250; 251 = reference self; 252 = not scored
  "islands":  [{"s":74,"e":109,"m":0.8233,"plp":3,"blb":0}],  // plp/blb = ClinVar residues inside
  "domains":  [{"name":"Kinesin motor","s":5,"e":354}],
  "clinvar":        [{"p":662,"ref":"P","alt":"L","sig":"P|PLP|LP","n":1}],  // up (pathogenic)
  "clinvar_benign": [{"p":..., "sig":"B|BLB|LB", ...}]                        // down (benign)
}
```

`index.json` = `[{"symbol","uniprot","len"}]`, drives search autocomplete.

## Precompute pipeline (`scripts/`, Python 3, run locally)

Run all: `python3 scripts/run_all.py [--genes KIF1A,TP53]` (omit `--genes` for full run).
Individual stages:

0. **`build_uniprot_map.py`** — one pass over `AlphaMissense_hg38.tsv.gz`; take each
   UniProt's genomic locus and overlap it with MANE gene spans → cache
   `scripts/.uniprot_map.json` = `{uniprot: {symbol, ensembl, chrom}}`. Needed because the
   aa_substitutions file has no symbol/coords, and AM's transcript is often **not** MANE
   Select (so an ID join misses ~28% of genes — overlap is required). Auto-run by run_all
   when the cache is absent.
1. **`build_bundles.py`** — stream `AlphaMissense_aa_substitutions.tsv.gz` grouped by
   UniProt (file is uniprot-sorted); per residue collect ref aa + all-19 scores; flush one
   bundle per protein. Symbol/ensembl from the UniProt map; islands + domains joined by
   UniProt. Emit per-gene bundles (no ClinVar yet).
2. **`add_clinvar.py`** — one pass over `clinvar.vep.vcf.gz`; take P/LP + B/LB missense,
   grouped by (symbol, RefSeq NM). Per gene **auto-select the RefSeq isoform** whose ref aa
   best matches the AM sequence (AM often uses a non-MANE isoform), collapse per residue to
   strongest call (P>PLP>LP; B>BLB>LB) with count `n`, drop off-isoform positions, set
   `refseq` to the chosen NM, and annotate each island's `plp`/`blb`.
3. **`build_index.py`** — emit `public/data/index.json`.

`scripts/config.py` holds source paths + constants (AA_ORDER, score cutoffs 0.34/0.564,
quantization). `scripts/lib.py` holds shared parsing (MANE map, base64 grid pack/unpack,
domain-range parser, CSQ parser).

Re-run when a source updates (new AM release, new ClinVar dump, new MANE). Bump the
provenance strings in the footer + this file.

## Frontend (static Vite SPA, no server)

- `src/main.ts` — landing↔viewer routing via `?gene=SYMBOL`; loads `index.json`; fetches +
  inflates the bundle; theme toggle (persisted).
- `src/search.ts` — reusable autocomplete widget (hero bar + mini top-bar bar).
- `src/viewer.ts` — canvas track stack (axis, domains, islands, ClinVar bidirectional
  lollipops, mean-AM line, continuous heatmap), zoom/fit, crosshair, tooltip.
- `src/styles.css` — design tokens + light/dark themes.

## Deploy (Vercel)

Framework preset: Vite. Build `npm run build` → `dist/`. `public/data/**` is copied
verbatim. `vercel.json` sets long cache on `/data/genes/**`. No env vars, no serverless.

## Provenance / caveats

- Heatmap = **all 19 substitutions** per residue (dense), from the aa_substitutions file.
  Continuous AM scores quantized to 1 byte (≈0.004 resolution). NB: the **islands** and the
  hg38-derived symbol map come from the SNV-based release — both are AlphaMissense, but the
  island means were computed over SNV-reachable substitutions, a slight source nuance vs the
  all-19 heatmap.
- ClinVar: **v107 (2025-10)**; only P/LP and B/LB missense; conflicting/VUS excluded
  (could be added as a toggle track later). Per-residue collapsed to the strongest call.
- ref-aa concordance gate protects against AM-isoform vs RefSeq mismatches.
- Repo carries ~19k gzipped bundles; if the plain commit gets heavy, move to Git LFS.
