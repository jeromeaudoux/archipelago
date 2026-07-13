# AM Islands Viewer

Type a human gene → explore its **AlphaMissense missense saturation map**, pathogenic
**islands** (colored by ClinVar overlap), **ClinVar** P/LP & B/LB variants, and **protein
domains** — all aligned on the protein sequence.

Static Vite SPA. Data is **precomputed offline** into one small gzipped JSON per gene and
served as static files, so the browser never touches the multi-hundred-MB source datasets.

## Develop

```bash
npm install
npm run dev            # http://localhost:5173  (open ?gene=KIF1A)
```

## Precompute the gene bundles

The raw datasets are large and live outside the repo — see paths in `scripts/config.py`.

```bash
python3 scripts/run_all.py --genes KIF1A,TP53   # quick spot-check
python3 scripts/run_all.py                       # full run (all MANE genes)
```

Outputs to `public/data/genes/<SYMBOL>.json.gz` and `public/data/index.json`.

## Build & deploy

```bash
npm run build          # → dist/
```

Deploy to Vercel (framework preset: Vite). `public/data/**` ships as static assets.

See `CLAUDE.md` for architecture, data sources, bundle schema, and provenance.
