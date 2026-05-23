#!/usr/bin/env python3
"""
00_make_manifest.py — Generate samples_100.tsv for the 100-sample chr22 test

Reads: config/samples.tsv (NHGRI_ID, Sex, SubPopulation, SuperPopulation, ONT_library, ONT_pore)
Does:  For each of the first N samples, lists S3 MINIMAP2_ALIGNED_BAMS/{sample_id}/ to
       find the actual BAM filename and derive bam_prefix.
Writes: config/samples_100.tsv in the same format as samples_pilot.tsv

Usage (run on Gadi login node — awscli works there):
    module load python3/3.12.1
    python3 scripts/00_pilot/00_make_manifest.py \\
        --samples-tsv  config/samples.tsv \\
        --out          /scratch/cy94/sb8857/1kgp-mqtl-100/config/samples_100.tsv \\
        --n-samples    100

Notes:
- bam_prefix = BAM filename without .phased.bam extension
- basecaller_model = 0 (R9) or 1 (R10) — binary encoding for tensorQTL covariate
- s3_bam_dir, s3_vcf_dir, s3_gvcf_dir are fixed paths (same for all samples)
- Requires: awscli installed (module load python3/3.12.1; pip install --user awscli)
"""

import argparse
import subprocess
import sys
import logging

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%H:%M:%S")
log = logging.getLogger(__name__)

S3_BASE = "s3://1000g-ont/PROCESSED_DATA/ALIGNED_TO_HG38"
# Full path to awscli — installed via 'pip install --user' under python3/3.12.1
# 'aws' alone fails in nohup/PBS shells where ~/.local/bin is not on PATH
import os as _os
AWS_CMD = _os.path.expanduser("~/.local/bin/aws")
S3_BAM_DIR = "MINIMAP2_ALIGNED_BAMS"
S3_VCF_DIR = "CLAIR3/PHASED_VCF"
S3_GVCF_DIR = "CLAIR3/GVCF"

OUT_HEADER = "sample_id\tsex\tsub_pop\tsuper_pop\tlibrary\tpore\tbasecaller_model\tbam_prefix\ts3_bam_dir\ts3_vcf_dir\ts3_gvcf_dir"


def get_bam_prefix(sample_id: str) -> str:
    """Find the .phased.bam for this sample and return its prefix.

    BAMs are stored at the root of MINIMAP2_ALIGNED_BAMS/ (NOT in per-sample
    subdirectories). E.g.:
      MINIMAP2_ALIGNED_BAMS/GM18507-ONT-hg38-R9-LSK110-dorado050_sup_5mCG_v33.phased.bam

    Strategy: prefix search on the sample_id (no trailing slash) so we only
    retrieve files belonging to this sample.

    Note: shell=True is required — Python subprocess.run with the aws binary
    fails (returncode=1) even when the shell command works, because the aws
    shebang (#!/usr/bin/env python3) doesn't resolve the module-loaded python3
    correctly when called via execve from within a subprocess. shell=True
    uses the already-configured bash environment. (Discovered 2026-03-14,
    copyq node debug job 163058967.)
    """
    # Prefix search — no trailing slash, no sample subdir
    s3_prefix = f"{S3_BASE}/{S3_BAM_DIR}/{sample_id}"
    cmd = f"{AWS_CMD} s3 ls {s3_prefix} --no-sign-request"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"S3 listing failed for {sample_id}: "
            f"stdout={result.stdout.strip()!r} stderr={result.stderr.strip()!r}"
        )

    for line in result.stdout.splitlines():
        parts = line.strip().split()
        if not parts:
            continue
        fname = parts[-1]
        if fname.endswith(".phased.bam") and not fname.endswith(".bai"):
            return fname.replace(".phased.bam", "")

    raise RuntimeError(
        f"No .phased.bam found for {sample_id} at prefix {s3_prefix}. "
        f"Raw output: {result.stdout.strip()!r}"
    )


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--samples-tsv", required=True, help="config/samples.tsv (500-sample master list)")
    p.add_argument("--out", required=True, help="Output manifest path")
    p.add_argument("--n-samples", type=int, default=100, help="Number of samples to include")
    return p.parse_args()


def main():
    args = parse_args()

    # Read master samples.tsv
    samples = []
    with open(args.samples_tsv) as f:
        header = f.readline().strip().split("\t")
        log.info(f"Columns: {header}")
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 6:
                continue
            samples.append({
                "sample_id": parts[0],
                "sex":        parts[1],
                "sub_pop":    parts[2],
                "super_pop":  parts[3],
                "library":    parts[4],
                "pore":       parts[5],
            })

    samples = samples[:args.n_samples]
    log.info(f"Processing {len(samples)} samples")

    rows = []
    failed = []
    for i, s in enumerate(samples):
        sid = s["sample_id"]
        try:
            bam_prefix = get_bam_prefix(sid)
            # basecaller_model: 0=R9, 1=R10 (binary covariate encoding, as in pilot)
            basecaller_model = 1 if s["pore"] == "R10" else 0
            row = (
                f"{sid}\t{s['sex']}\t{s['sub_pop']}\t{s['super_pop']}\t"
                f"{s['library']}\t{s['pore']}\t{basecaller_model}\t"
                f"{bam_prefix}\t{S3_BAM_DIR}\t{S3_VCF_DIR}\t{S3_GVCF_DIR}"
            )
            rows.append(row)
            log.info(f"[{i+1}/{len(samples)}] {sid} → {bam_prefix}")
        except Exception as e:
            log.error(f"[{i+1}/{len(samples)}] {sid} FAILED: {e}")
            failed.append(sid)

    if failed:
        log.error(f"Failed samples ({len(failed)}): {failed}")
        log.error("Fix failures before proceeding — manifest must cover all samples")
        sys.exit(1)

    # Write output
    import os
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as f:
        f.write(OUT_HEADER + "\n")
        for row in rows:
            f.write(row + "\n")

    log.info(f"Written: {args.out} ({len(rows)} samples)")

    # Summary: pore breakdown (verify R9/R10 split)
    r9 = sum(1 for r in rows if "\tR9\t" in r)
    r10 = sum(1 for r in rows if "\tR10\t" in r)
    log.info(f"Pore breakdown: R9={r9}, R10={r10}")


if __name__ == "__main__":
    main()
