"""Stage 2 — merge ClinVar P/LP + B/LB missense into the gene bundles.

One streaming pass over clinvar.vep.vcf.gz collects, per (gene symbol, RefSeq NM),
the ClinVar missense calls. For each gene we then auto-select the RefSeq transcript
whose reference amino acids are most concordant with the AlphaMissense isoform (AM
often uses a non-MANE isoform), attach that transcript's variants, and annotate each
island's P/LP and B/LB overlap counts.

Usage:
  python3 add_clinvar.py [--genes KIF1A,TP53]
"""
import argparse
import glob
import gzip
import json
import os
import re
import sys
import time
from collections import defaultdict

import config as C
import lib

CLNSIG_RE = re.compile(r"CLNSIG=([^;]+)")
CSQ_RE = re.compile(r"CSQ=([^;\t]+)")

LOF_TERMS = ("stop_gained", "frameshift_variant", "splice_acceptor_variant",
             "splice_donor_variant", "start_lost", "stop_lost", "transcript_ablation")


def categorize(cons):
    if "missense" in cons:
        return "missense"
    if any(t in cons for t in LOF_TERMS):
        return "lof"
    return "other"


def parse_csq_format(vcf_path):
    with lib.openf(vcf_path) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            if line.startswith("##INFO=<ID=CSQ"):
                fmt = line.split("Format:", 1)[1].strip().rstrip('">').strip()
                return {name: i for i, name in enumerate(fmt.split("|"))}
    raise RuntimeError("CSQ format not found in VCF header")


def collect_clinvar(vcf_path):
    ci = parse_csq_format(vcf_path)
    i_sym, i_feat = ci["SYMBOL"], ci["Feature"]
    i_cons, i_ppos, i_aa = ci["Consequence"], ci["Protein_position"], ci["Amino_acids"]
    ncols = max(i_sym, i_feat, i_cons, i_ppos, i_aa) + 1

    # symbol -> nm_base -> pos -> rec  (missense only; drives the heatmap track)
    patho = defaultdict(lambda: defaultdict(dict))
    benign = defaultdict(lambda: defaultdict(dict))
    nm_full = {}  # (symbol, nm_base) -> versioned NM
    # (symbol, nm_base) -> {missense|lof|other: [P/LP count, B/LB count]} — all consequences
    breakdown = defaultdict(lambda: {"missense": [0, 0], "lof": [0, 0], "other": [0, 0]})
    t0 = time.time()
    seen = 0

    with lib.openf(vcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            m = CLNSIG_RE.search(line)
            if not m:
                continue
            clnsig = m.group(1)
            if clnsig in C.PATHO_SIG:
                sig, rank, store, pb = C.PATHO_SIG[clnsig], C.PATHO_RANK, patho, 0
            elif clnsig in C.BENIGN_SIG:
                sig, rank, store, pb = C.BENIGN_SIG[clnsig], C.BENIGN_RANK, benign, 1
            else:
                continue
            mc = CSQ_RE.search(line)
            if not mc:
                continue
            seen += 1
            for blk in mc.group(1).split(","):
                p = blk.split("|")
                if len(p) < ncols:
                    continue
                symbol, feat, cons = p[i_sym], p[i_feat], p[i_cons]
                if not symbol or not feat.startswith("NM_"):
                    continue
                nm_base = lib.base(feat)
                breakdown[(symbol, nm_base)][categorize(cons)][pb] += 1
                if "missense" not in cons or "/" not in p[i_aa]:
                    continue
                try:
                    pos = int(p[i_ppos].split("-")[0])
                except ValueError:
                    continue
                nm_full.setdefault((symbol, nm_base), feat)
                aa = p[i_aa]
                ref, alt = aa.split("/")[0], aa.split("/")[-1]
                d = store[symbol][nm_base]
                cur = d.get(pos)
                if cur is None:
                    d[pos] = {"p": pos, "ref": ref, "alt": alt, "sig": sig, "n": 1}
                else:
                    cur["n"] += 1
                    if rank[sig] > rank[cur["sig"]]:
                        cur["sig"], cur["alt"] = sig, alt

    print(f"  ClinVar scan: {seen} P/LP+B/LB records ({time.time()-t0:.0f}s)", flush=True)
    return patho, benign, nm_full, breakdown


def concordance(nm_variants, ref_str, length):
    ok = 0
    for pos, rec in nm_variants.items():
        if 1 <= pos <= length and ref_str[pos - 1] == rec["ref"]:
            ok += 1
    return ok


def resolve(nm_variants, ref_str, length):
    out, mism = [], 0
    for pos, rec in sorted(nm_variants.items()):
        if 1 <= pos <= length and ref_str[pos - 1] == rec["ref"]:
            out.append(rec)
        else:
            mism += 1
    return out, mism


def inject(patho, benign, nm_full, breakdown, targets):
    files = sorted(glob.glob(os.path.join(C.GENES_DIR, "*.json.gz")))
    updated = mism_total = 0
    for path in files:
        symbol = os.path.basename(path)[: -len(".json.gz")]
        if targets is not None and symbol not in targets:
            continue
        with gzip.open(path, "rt") as fh:
            b = json.load(fh)
        ref_str, length = b["ref"], b["length"]

        # Candidate RefSeq NMs for this gene from either bucket; pick the isoform
        # whose ref aa best matches the AlphaMissense sequence.
        candidates = set(patho.get(symbol, {})) | set(benign.get(symbol, {}))
        best_nm, best_score = None, -1
        for nm in candidates:
            score = (concordance(patho.get(symbol, {}).get(nm, {}), ref_str, length)
                     + concordance(benign.get(symbol, {}).get(nm, {}), ref_str, length))
            if score > best_score:
                best_nm, best_score = nm, score

        cv, cvb = [], []
        if best_nm is not None:
            cv, m1 = resolve(patho.get(symbol, {}).get(best_nm, {}), ref_str, length)
            cvb, m2 = resolve(benign.get(symbol, {}).get(best_nm, {}), ref_str, length)
            mism_total += m1 + m2
            b["refseq"] = nm_full.get((symbol, best_nm), best_nm)
        b["clinvar"], b["clinvar_benign"] = cv, cvb

        # Consequence breakdown (missense/LoF/other) on the chosen transcript; if the
        # gene has no missense ClinVar, fall back to its most-populated transcript.
        bd_nm = best_nm
        if bd_nm is None:
            cands = {nm for (s, nm) in breakdown if s == symbol}
            if cands:
                bd_nm = max(cands, key=lambda nm: sum(
                    sum(v) for v in breakdown[(symbol, nm)].values()))
        bd = breakdown.get((symbol, bd_nm)) if bd_nm is not None else None
        if bd and any(sum(v) for v in bd.values()):
            b["cv_counts"] = bd

        ppos = {r["p"] for r in cv}
        bpos = {r["p"] for r in cvb}
        for isl in b["islands"]:
            isl["plp"] = sum(1 for p in ppos if isl["s"] <= p <= isl["e"])
            isl["blb"] = sum(1 for p in bpos if isl["s"] <= p <= isl["e"])

        payload = json.dumps(b, separators=(",", ":")).encode("utf-8")
        with gzip.GzipFile(path, "wb", mtime=0) as gz:
            gz.write(payload)
        updated += 1
    print(f"  Injected ClinVar into {updated} bundles "
          f"({mism_total} off-isoform positions dropped)", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", help="comma-separated gene symbols to limit")
    args = ap.parse_args()
    targets = set(args.genes.split(",")) if args.genes else None

    if not os.path.isdir(C.GENES_DIR):
        print("No bundles found — run build_bundles.py first.", file=sys.stderr)
        sys.exit(1)
    print(f"Scanning {C.CLINVAR_VCF}…", flush=True)
    patho, benign, nm_full, breakdown = collect_clinvar(C.CLINVAR_VCF)
    inject(patho, benign, nm_full, breakdown, targets)


if __name__ == "__main__":
    main()
