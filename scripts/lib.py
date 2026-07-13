"""Shared parsing helpers for the precompute pipeline."""
import base64
import gzip
import re
from collections import defaultdict

import config as C


def openf(path):
    """Open a text file, transparently gunzipping .gz."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def base(feature_id):
    """Strip a version suffix: ENST00000650053.1 -> ENST00000650053."""
    return feature_id.split(".", 1)[0] if feature_id else feature_id


# ---------------------------------------------------------------------------
# MANE Select mapping backbone
# ---------------------------------------------------------------------------
def load_mane(select_only=True):
    """Parse the MANE summary.

    Returns (by_enst, by_symbol) where each maps to a dict with:
      symbol, enst_base, enst_full, refseq_base, refseq_full.
    Only MANE Select rows are kept by default (one transcript per gene).
    """
    by_enst, by_symbol = {}, {}
    with openf(C.MANE_SUMMARY) as fh:
        header = fh.readline().rstrip("\n").lstrip("#").split("\t")
        col = {name: i for i, name in enumerate(header)}
        for line in fh:
            v = line.rstrip("\n").split("\t")
            status = v[col["MANE_status"]]
            if select_only and status != "MANE Select":
                continue
            enst_full = v[col["Ensembl_nuc"]]
            refseq_full = v[col["RefSeq_nuc"]]
            rec = {
                "symbol": v[col["symbol"]],
                "enst_base": base(enst_full),
                "enst_full": enst_full,
                "refseq_base": base(refseq_full),
                "refseq_full": refseq_full,
            }
            by_enst[rec["enst_base"]] = rec
            by_symbol[rec["symbol"]] = rec
    return by_enst, by_symbol


# ---------------------------------------------------------------------------
# AM islands (keyed by Ensembl transcript)
# ---------------------------------------------------------------------------
def _load_islands(key_col):
    """key_col: 3 = transcript_id (ENST), 4 = uniprot_id."""
    islands = defaultdict(list)
    with openf(C.AM_ISLANDS_BED) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            k = base(f[key_col]) if key_col == 3 else f[key_col]
            islands[k].append(
                {"s": int(f[7]), "e": int(f[8]), "m": round(float(f[5]), 4)}
            )
    for v in islands.values():
        v.sort(key=lambda d: d["s"])
    return islands


def load_islands():
    """enst_base -> list of {s, e, m}."""
    return _load_islands(3)


def load_islands_by_uniprot():
    """uniprot -> list of {s, e, m}."""
    return _load_islands(4)


# ---------------------------------------------------------------------------
# UniProt domains (keyed by UniProt accession)
# ---------------------------------------------------------------------------
_DOMAIN_RE = re.compile(r"^([A-Z])(\d+)-([A-Z])(\d+);\s*(.+)$")


def load_domains():
    """uniprot -> list of {name, s, e} (deduped, sorted)."""
    domains = defaultdict(list)
    seen = defaultdict(set)
    with openf(C.DOMAIN_BED) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 14:
                continue
            uniprot = f[3]
            m = _DOMAIN_RE.match(f[13].strip())
            if not m:
                continue
            s, e, name = int(m.group(2)), int(m.group(4)), m.group(5).strip()
            key = (s, e, name)
            if key in seen[uniprot]:
                continue
            seen[uniprot].add(key)
            domains[uniprot].append({"name": name, "s": s, "e": e})
    for v in domains.values():
        v.sort(key=lambda d: d["s"])
    return domains


# ---------------------------------------------------------------------------
# Continuous-score grid packing
# ---------------------------------------------------------------------------
def pack_grid(length, ref_by_pos, scores_by_pos):
    """Build the base64 N*20 byte grid.

    ref_by_pos:    {pos(1-based) -> ref aa}
    scores_by_pos: {pos(1-based) -> {alt aa -> score float}}
    """
    buf = bytearray(length * len(C.AA_ORDER))
    for pos in range(1, length + 1):
        ref = ref_by_pos.get(pos)
        scores = scores_by_pos.get(pos, {})
        row = (pos - 1) * len(C.AA_ORDER)
        for aa, idx in C.AA_INDEX.items():
            if aa == ref:
                buf[row + idx] = C.BYTE_REFERENCE
            elif aa in scores:
                b = int(round(scores[aa] * C.Q_MAX))
                buf[row + idx] = 0 if b < 0 else (C.Q_MAX if b > C.Q_MAX else b)
            else:
                buf[row + idx] = C.BYTE_MISSING
    return base64.b64encode(bytes(buf)).decode("ascii")


# ---------------------------------------------------------------------------
# protein_variant parsing:  "V2L" -> ("V", 2, "L")
# ---------------------------------------------------------------------------
_PVAR_RE = re.compile(r"^([A-Z\*])(\d+)([A-Z\*])$")


def parse_pvar(pvar):
    m = _PVAR_RE.match(pvar)
    if not m:
        return None
    return m.group(1), int(m.group(2)), m.group(3)


# ---------------------------------------------------------------------------
# MANE gene intervals — map an AlphaMissense transcript to a gene symbol by
# genomic overlap (AM's transcript is often NOT the MANE Select one, so an
# ID join misses ~28% of genes; a coordinate overlap is robust).
# ---------------------------------------------------------------------------
def _nc_to_chr(acc):
    """NC_000002.12 -> chr2 ; NC_000023 -> chrX ; NC_012920 -> chrM."""
    try:
        num = int(acc.split(".")[0].split("_")[1])
    except (IndexError, ValueError):
        return None
    if 1 <= num <= 22:
        return "chr" + str(num)
    if num == 23:
        return "chrX"
    if num == 24:
        return "chrY"
    if num == 12920:
        return "chrM"
    return None


def load_mane_intervals(select_only=True):
    """chrom -> list of (start, end, symbol), sorted by start."""
    per = defaultdict(list)
    with openf(C.MANE_SUMMARY) as fh:
        header = fh.readline().rstrip("\n").lstrip("#").split("\t")
        col = {n: i for i, n in enumerate(header)}
        for line in fh:
            v = line.rstrip("\n").split("\t")
            if select_only and v[col["MANE_status"]] != "MANE Select":
                continue
            chrom = _nc_to_chr(v[col["GRCh38_chr"]])
            if not chrom:
                continue
            per[chrom].append(
                (int(v[col["chr_start"]]), int(v[col["chr_end"]]), v[col["symbol"]])
            )
    for v in per.values():
        v.sort(key=lambda t: t[0])
    return per


def lookup_symbol(intervals, chrom, pos):
    """Symbol of the smallest MANE gene interval containing pos (or None)."""
    lst = intervals.get(chrom)
    if not lst:
        return None
    best, best_span = None, None
    for s, e, sym in lst:
        if s <= pos <= e:
            span = e - s
            if best_span is None or span < best_span:
                best, best_span = sym, span
        elif s > pos:
            break
    return best
