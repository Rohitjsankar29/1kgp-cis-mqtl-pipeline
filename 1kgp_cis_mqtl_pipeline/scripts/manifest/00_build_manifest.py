#!/usr/bin/env python3
"""
Build a 1KGP ONT BAM manifest for a small pilot or full run.

Output columns:
sample_id, bam_s3_url, bai_s3_url, bam_filename, bam_size_gb,
pore, basecaller_model, basecaller, modifications, has_5hmc,
sex, sub_pop, super_pop, library

This avoids downloading BAMs during manifest creation.
"""
import argparse, csv, re, sys, urllib.parse, urllib.request, xml.etree.ElementTree as ET
from pathlib import Path

S3_HTTP = "https://s3.amazonaws.com/1000g-ont"
BAM_PREFIX = "PROCESSED_DATA/ALIGNED_TO_HG38/MINIMAP2_ALIGNED_BAMS/"
META_URL = "https://s3.amazonaws.com/1000g-ont/PROCESSED_DATA/1kGP_LRSC_500_ONT_Metadata.tsv"
NS = "http://s3.amazonaws.com/doc/2006-03-01/"

def list_s3_prefix(prefix):
    out, token = [], None
    while True:
        url = f"{S3_HTTP}?list-type=2&prefix={urllib.parse.quote(prefix)}&max-keys=1000"
        if token:
            url += f"&continuation-token={urllib.parse.quote(token)}"
        with urllib.request.urlopen(url) as r:
            root = ET.parse(r).getroot()
        for c in root.findall(f"{{{NS}}}Contents"):
            key = c.find(f"{{{NS}}}Key").text
            size = int(c.find(f"{{{NS}}}Size").text)
            out.append((key, size))
        trunc = root.find(f"{{{NS}}}IsTruncated")
        if trunc is not None and trunc.text == "true":
            token = root.find(f"{{{NS}}}NextContinuationToken").text
        else:
            return out

def parse_bam(fname):
    m = re.match(r"^([A-Z]+\d+)-ONT-", fname)
    sample_id = m.group(1) if m else fname.split("-ONT-")[0]
    pore = "R10" if "-R10-" in fname else "R9"
    dm = re.search(r"dorado(\d+)", fname)
    if "guppy657" in fname: basecaller = "guppy657"
    elif "guppy" in fname: basecaller = "guppy"
    elif dm: basecaller = "dorado" + dm.group(1)
    else: basecaller = "unknown"
    mods = ["5mC"]
    if "5hmCG" in fname or "5hmC" in fname:
        mods.append("5hmC")
    return sample_id, pore, basecaller, "+".join(mods), ("5hmC" in mods)

def load_metadata(path_or_url):
    meta = {}
    try:
        fh = urllib.request.urlopen(path_or_url) if str(path_or_url).startswith("http") else open(path_or_url, "rb")
        text = fh.read().decode().splitlines()
        header = text[0].split("\t")
        lower = [h.lower() for h in header]
        for line in text[1:]:
            vals = line.split("\t")
            row = dict(zip(lower, vals))
            sid = row.get("nhgri_id") or row.get("sample_id") or vals[0]
            meta[sid] = row
    except Exception as e:
        print(f"WARN: metadata unavailable: {e}", file=sys.stderr)
    return meta

def pick(row, names, default=""):
    for n in names:
        if n in row and row[n] != "": return row[n]
    return default

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="config/bam_manifest.tsv")
    ap.add_argument("--n-samples", type=int, default=2)
    ap.add_argument("--samples", default="", help="Optional comma-separated sample IDs to keep")
    ap.add_argument("--metadata", default=META_URL)
    args = ap.parse_args()

    wanted = set(x.strip() for x in args.samples.split(",") if x.strip())
    meta = load_metadata(args.metadata)
    objects = list_s3_prefix(BAM_PREFIX)
    bams = [(k,s) for k,s in objects if k.endswith(".phased.bam")]
    rows = []
    for key, size in sorted(bams):
        fname = Path(key).name
        sid, pore, basecaller, mods, has_5hmc = parse_bam(fname)
        if wanted and sid not in wanted: continue
        m = meta.get(sid, {})
        rows.append({
            "sample_id": sid,
            "bam_s3_url": f"s3://1000g-ont/{key}",
            "bai_s3_url": f"s3://1000g-ont/{key}.bai",
            "bam_filename": fname,
            "bam_size_gb": round(size / 1e9, 2),
            "pore": pore,
            "basecaller_model": 1 if pore == "R10" else 0,
            "basecaller": basecaller,
            "modifications": mods,
            "has_5hmc": str(has_5hmc).upper(),
            "sex": pick(m, ["sex", "gender"], "NA"),
            "sub_pop": pick(m, ["subpopulation", "sub_pop", "population"], "NA"),
            "super_pop": pick(m, ["superpopulation", "super_pop"], "NA"),
            "library": pick(m, ["ont_library", "library"], "NA"),
        })
        if not wanted and len(rows) >= args.n_samples: break
    if not rows: sys.exit("No BAM rows found")
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(rows)
    print(f"Wrote {len(rows)} samples: {args.out}", file=sys.stderr)

if __name__ == "__main__": main()
