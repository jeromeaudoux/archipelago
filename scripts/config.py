"""Shared configuration for the precompute pipeline.

Raw source datasets are large and live OUTSIDE the repo. Adjust paths here if your
local layout differs. Nothing in here is shipped to the browser.
"""
import os

HOME = os.path.expanduser("~")

# --- Raw inputs (not committed) -------------------------------------------------
ALPHAMISSENSE_TSV = os.path.join(HOME, "Downloads", "AlphaMissense_hg38.tsv.gz")
# All-19-substitutions release (keyed by UniProt; dense heatmap source).
ALPHAMISSENSE_AASUB = os.path.join(HOME, "Downloads", "AlphaMissense_aa_substitutions.tsv.gz")
# Cache: uniprot -> {symbol, ensembl, chrom} derived from the hg38 file (genomic
# overlap with MANE gene spans), so the UniProt-keyed aa_substitutions data can be
# named by gene symbol.
UNIPROT_MAP_JSON = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".uniprot_map.json")
AM_ISLANDS_BED = os.path.join(HOME, "Dev", "acmgscore-v2", "data", "am_islands_hg38.bed.gz")
DOMAIN_BED = os.path.join(HOME, "Dev", "acmgscore-v2", "UP000005640_9606_domain.bed")
CLINVAR_VCF = os.path.join(
    HOME, "Dev", "Data", "resources", "human", "clinvar-vep", "GRCh38",
    "107_2025-10", "clinvar.vep.vcf.gz",
)
MANE_SUMMARY = os.path.join(HOME, "Dev", "Data", "MANE.GRCh38.v1.3.summary.txt.gz")
# cancerhotspots.org residue hotspots (hotspots_v3.xlsx, "Hotspot_Residues" sheet).
CANCERHOTSPOTS_XLSX = os.path.join(HOME, "Downloads", "hotspots_v3.xlsx")
CLINVAR_VERSION = "v107 · 2025-10"

# --- Outputs --------------------------------------------------------------------
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(REPO, "public", "data")
GENES_DIR = os.path.join(DATA_DIR, "genes")
INDEX_JSON = os.path.join(DATA_DIR, "index.json")

# --- Model constants ------------------------------------------------------------
# Heatmap row order: charged+ / charged- / polar / special / hydrophobic
AA_ORDER = "HKRDESTNQCGPAVILMFYW"
AA_INDEX = {a: i for i, a in enumerate(AA_ORDER)}

# AlphaMissense class boundaries (for reference / UI guides).
CUTOFF_BENIGN = 0.34
CUTOFF_PATHO = 0.564

# Continuous-score quantization for the packed grid.
Q_MAX = 250          # score 0..1 -> byte 0..250
BYTE_REFERENCE = 251  # the reference amino acid at that position
BYTE_MISSING = 252    # substitution not scored by AlphaMissense

# ClinVar significance buckets.
PATHO_SIG = {
    "Pathogenic": "P",
    "Pathogenic/Likely_pathogenic": "PLP",
    "Likely_pathogenic": "LP",
}
BENIGN_SIG = {
    "Benign": "B",
    "Benign/Likely_benign": "BLB",
    "Likely_benign": "LB",
}
PATHO_RANK = {"P": 3, "PLP": 2, "LP": 1}
BENIGN_RANK = {"B": 3, "BLB": 2, "LB": 1}
