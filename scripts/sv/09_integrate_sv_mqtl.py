#!/usr/bin/env python3
"""
09_integrate_sv_mqtl.py

Integrate structural variants into the CpG-centric mQTL feature table.

Purpose:
    The dissertation framework includes SNPs and structural variants (SVs) as
    candidate regulatory variants. This script annotates variant-CpG pairs with
    SV evidence and optionally creates SV-CpG cis candidate pairs.

Inputs:
    --features      Existing feature matrix TSV from SNP cis-mQTL results
    --sv-vcf        SV VCF/VCF.GZ from ONT/Sniffles or other caller
    --cpg-bed       CpG methylation phenotype BED.gz (#chr start end phenotype_id ...)
    --out           Updated feature matrix TSV

Outputs:
    Adds columns:
        is_sv
        sv_type
        sv_len
        sv_overlaps_cpg
        sv_distance_to_cpg
        sv_within_cis_window

Notes:
    - This script does not perform SV association testing itself.
    - It creates annotation features so SV evidence can enter prioritisation.
    - If you have SV-mQTL results, merge them first or pass them as --sv-qtl.

Author: Rohit pipeline helper
"""

import argparse
import gzip
import logging
import re
from pathlib import Path

import pandas as pd
import numpy as np

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%H:%M:%S")
log = logging.getLogger(__name__)


def open_text(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def parse_info(info):
    d = {}
    for item in info.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
        else:
            d[item] = True
    return d


def load_sv_vcf(vcf_path):
    rows = []
    with open_text(vcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, vid, ref, alt, qual, filt, info = parts[:8]
            pos = int(pos)
            info_d = parse_info(info)
            svtype = info_d.get("SVTYPE", "SV")
            end = int(info_d.get("END", pos))
            svlen = info_d.get("SVLEN", end - pos + 1)
            try:
                svlen = int(str(svlen).split(",")[0])
            except Exception:
                svlen = end - pos + 1
            if vid == ".":
                vid = f"{chrom}:{pos}:{svtype}:{abs(svlen)}"
            rows.append({
                "variant_id": vid,
                "chr": chrom if str(chrom).startswith("chr") else f"chr{chrom}",
                "sv_start": pos - 1,
                "sv_end": end,
                "sv_type": svtype,
                "sv_len": svlen,
            })
    if not rows:
        raise ValueError(f"No SV records found in {vcf_path}")
    sv = pd.DataFrame(rows)
    log.info(f"Loaded SVs: {len(sv):,}")
    return sv


def load_cpg_bed(cpg_bed):
    df = pd.read_csv(cpg_bed, sep="\t", compression="infer")
    chrom_col = "#chr" if "#chr" in df.columns else "chr"
    required = {chrom_col, "start", "end", "phenotype_id"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"CpG BED missing columns: {missing}")
    out = df[[chrom_col, "start", "end", "phenotype_id"]].copy()
    out = out.rename(columns={chrom_col: "chr", "start": "cpg_start", "end": "cpg_end"})
    out["cpg_mid"] = ((out["cpg_start"] + out["cpg_end"]) / 2).astype(int)
    out["chr"] = out["chr"].astype(str)
    log.info(f"Loaded CpGs: {len(out):,}")
    return out


def parse_variant_position(variant_id):
    # Supports chr:pos:ref:alt, chr_pos_ref_alt, chr:pos
    s = str(variant_id)
    m = re.match(r"^(chr)?([0-9XYM]+)[:_](\d+)", s)
    if not m:
        return (np.nan, np.nan)
    chrom = "chr" + m.group(2) if not str(m.group(2)).startswith("chr") else m.group(2)
    return chrom, int(m.group(3))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--features", required=True, help="Existing feature matrix TSV")
    p.add_argument("--sv-vcf", required=True, help="SV VCF/VCF.GZ")
    p.add_argument("--cpg-bed", required=True, help="Methylation phenotype BED.gz")
    p.add_argument("--out", required=True)
    p.add_argument("--cis-window", type=int, default=1_000_000)
    p.add_argument("--emit-sv-candidates", action="store_true",
                   help="Also append SV-CpG candidate pairs within cis window")
    args = p.parse_args()

    features = pd.read_csv(args.features, sep="\t")
    sv = load_sv_vcf(args.sv_vcf)
    cpg = load_cpg_bed(args.cpg_bed)

    # Default SV columns for existing SNP rows
    features["is_sv"] = features.get("is_sv", 0)
    features["sv_type"] = features.get("sv_type", "")
    features["sv_len"] = features.get("sv_len", np.nan)
    features["sv_overlaps_cpg"] = features.get("sv_overlaps_cpg", 0)
    features["sv_distance_to_cpg"] = features.get("sv_distance_to_cpg", np.nan)
    features["sv_within_cis_window"] = features.get("sv_within_cis_window", 0)

    # If existing rows contain SV variant IDs, annotate them
    sv_small = sv[["variant_id", "sv_type", "sv_len"]].drop_duplicates("variant_id")
    features = features.merge(sv_small.assign(is_sv_match=1), on="variant_id", how="left")
    m = features["is_sv_match"].fillna(0).astype(int) == 1
    features.loc[m, "is_sv"] = 1
    features.loc[m, "sv_type"] = features.loc[m, "sv_type_y"] if "sv_type_y" in features.columns else features.loc[m, "sv_type"]
    features.loc[m, "sv_len"] = features.loc[m, "sv_len_y"] if "sv_len_y" in features.columns else features.loc[m, "sv_len"]

    for col in ["sv_type_x", "sv_type_y", "sv_len_x", "sv_len_y", "is_sv_match"]:
        if col in features.columns:
            features = features.drop(columns=[col])

    if args.emit_sv_candidates:
        log.info("Creating SV-CpG candidate pairs within cis window")
        candidate_rows = []
        for chrom, cpg_chr in cpg.groupby("chr"):
            sv_chr = sv[sv["chr"] == chrom]
            if sv_chr.empty:
                continue
            for _, cg in cpg_chr.iterrows():
                start = cg["cpg_mid"] - args.cis_window
                end = cg["cpg_mid"] + args.cis_window
                nearby = sv_chr[(sv_chr["sv_end"] >= start) & (sv_chr["sv_start"] <= end)].copy()
                if nearby.empty:
                    continue
                nearby["phenotype_id"] = cg["phenotype_id"]
                nearby["cpg_chr"] = chrom
                nearby["cpg_start"] = cg["cpg_start"]
                nearby["cpg_end"] = cg["cpg_end"]
                nearby["distance_to_cpg"] = np.minimum(
                    abs(nearby["sv_start"] - cg["cpg_mid"]),
                    abs(nearby["sv_end"] - cg["cpg_mid"])
                )
                nearby["is_sv"] = 1
                nearby["sv_overlaps_cpg"] = ((nearby["sv_start"] <= cg["cpg_end"]) & (nearby["sv_end"] >= cg["cpg_start"])).astype(int)
                nearby["sv_distance_to_cpg"] = nearby["distance_to_cpg"]
                nearby["sv_within_cis_window"] = 1
                candidate_rows.append(nearby[[
                    "phenotype_id", "variant_id", "distance_to_cpg",
                    "is_sv", "sv_type", "sv_len", "sv_overlaps_cpg",
                    "sv_distance_to_cpg", "sv_within_cis_window"
                ]])
        if candidate_rows:
            sv_candidates = pd.concat(candidate_rows, ignore_index=True)
            # Add missing feature columns for compatibility
            for col in features.columns:
                if col not in sv_candidates.columns:
                    sv_candidates[col] = np.nan
            sv_candidates = sv_candidates[features.columns]
            features = pd.concat([features, sv_candidates], ignore_index=True)
            log.info(f"Appended SV candidate rows: {len(sv_candidates):,}")

    features.to_csv(args.out, sep="\t", index=False)
    log.info(f"Written: {args.out} ({len(features):,} rows)")


if __name__ == "__main__":
    main()
