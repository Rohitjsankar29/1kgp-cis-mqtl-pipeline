#!/usr/bin/env python3
"""
make_manifest_parallel.py — parallel manifest generator for 500-sample run.
Re-uses existing resolved entries from samples_100.tsv, parallel S3 lookups for the rest.
"""
import argparse, subprocess, sys, logging, os
from concurrent.futures import ThreadPoolExecutor, as_completed

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%H:%M:%S")
log = logging.getLogger(__name__)

S3_BASE = "s3://1000g-ont/PROCESSED_DATA/ALIGNED_TO_HG38"
AWS_CMD = os.path.expanduser("~/.local/bin/aws")

OUT_HEADER = "sample_id\tsex\tsub_pop\tsuper_pop\tlibrary\tpore\tbasecaller_model\tbam_prefix\ts3_bam_dir\ts3_vcf_dir\ts3_gvcf_dir"

def get_bam_prefix(sample_id):
    s3_prefix = f"{S3_BASE}/MINIMAP2_ALIGNED_BAMS/{sample_id}"
    cmd = f"{AWS_CMD} s3 ls {s3_prefix} --no-sign-request"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"S3 failed for {sample_id}: {result.stderr.strip()!r}")
    for line in result.stdout.splitlines():
        parts = line.strip().split()
        if parts and parts[-1].endswith(".phased.bam"):
            return parts[-1].replace(".phased.bam", "")
    raise RuntimeError(f"No .phased.bam found for {sample_id}")

def parse_existing(path):
    known = {}
    if not os.path.exists(path):
        return known
    with open(path) as f:
        next(f)  # skip header
        for line in f:
            cols = line.rstrip('\n').split('\t')
            if len(cols) >= 8:
                known[cols[0]] = cols  # key by sample_id
    return known

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--samples-tsv", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--n-samples", type=int, default=500)
    p.add_argument("--existing", default="")
    p.add_argument("--workers", type=int, default=20)
    args = p.parse_args()

    # Load existing resolved entries
    known = {}
    if args.existing:
        known = parse_existing(args.existing)
        log.info(f"Loaded {len(known)} existing entries from {args.existing}")

    # Load input samples
    with open(args.samples_tsv) as f:
        header = f.readline()
        samples = [line.rstrip('\n').split('\t') for line in f]
    samples = samples[:args.n_samples]
    log.info(f"Processing {len(samples)} samples")

    # Determine which need S3 lookup
    to_resolve = [(i, s) for i, s in enumerate(samples) if s[0] not in known]
    log.info(f"Need S3 lookup: {len(to_resolve)} | Already known: {len(known)}")

    # Parallel S3 lookups
    results = {}
    errors = []
    def lookup(idx_sample):
        idx, s = idx_sample
        sample_id = s[0]
        try:
            prefix = get_bam_prefix(sample_id)
            return idx, sample_id, prefix, None
        except Exception as e:
            return idx, sample_id, None, str(e)

    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(lookup, t): t for t in to_resolve}
        done = 0
        for fut in as_completed(futs):
            idx, sid, prefix, err = fut.result()
            done += 1
            if err:
                log.error(f"[{done}/{len(to_resolve)}] {sid}: {err}")
                errors.append(sid)
            else:
                results[sid] = (idx, prefix)
                log.info(f"[{done}/{len(to_resolve)}] {sid} → {prefix}")

    if errors:
        log.error(f"{len(errors)} samples failed S3 lookup: {errors}")
        sys.exit(1)

    # Write output
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, 'w') as fout:
        fout.write(OUT_HEADER + '\n')
        for s in samples:
            sample_id = s[0]
            sex = s[1] if len(s) > 1 else ''
            sub_pop = s[2] if len(s) > 2 else ''
            super_pop = s[3] if len(s) > 3 else ''
            library = s[4] if len(s) > 4 else ''
            pore = s[5] if len(s) > 5 else ''

            if sample_id in known:
                fout.write('\t'.join(known[sample_id]) + '\n')
            else:
                prefix = results[sample_id][1]
                # basecaller_model: 0=R9, 1=R10
                basecaller_model = '1' if pore == 'R10' else '0'
                row = [sample_id, sex, sub_pop, super_pop, library, pore,
                       basecaller_model, prefix,
                       'MINIMAP2_ALIGNED_BAMS', 'CLAIR3/PHASED_VCF', 'CLAIR3/GVCF']
                fout.write('\t'.join(row) + '\n')

    log.info(f"Written {len(samples)} samples to {args.out}")

if __name__ == '__main__':
    main()
