#!/usr/bin/env python3
"""
Build bam_manifest.tsv: sample_id \t s3_bam_url \t basecaller \t modifications \t pore

Queries the S3 bucket, matches BAMs to metadata, annotates with
basecaller and modification model from the filename.

Output: config/bam_manifest.tsv
"""

import urllib.request
import xml.etree.ElementTree as ET
import csv
import re
import sys
from pathlib import Path

S3_BASE = "https://s3.amazonaws.com/1000g-ont"
BAM_PREFIX = "PROCESSED_DATA/ALIGNED_TO_HG38/MINIMAP2_ALIGNED_BAMS/"

NS = "http://s3.amazonaws.com/doc/2006-03-01/"


def list_s3_prefix(prefix: str) -> list[dict]:
    """List all objects under a given S3 prefix."""
    objects = []
    continuation = None
    while True:
        url = f"{S3_BASE}?list-type=2&prefix={urllib.parse.quote(prefix)}&max-keys=1000"
        if continuation:
            url += f"&continuation-token={urllib.parse.quote(continuation)}"
        with urllib.request.urlopen(url) as r:
            root = ET.parse(r).getroot()
        for c in root.findall(f"{{{NS}}}Contents"):
            key = c.find(f"{{{NS}}}Key").text
            size = int(c.find(f"{{{NS}}}Size").text)
            objects.append({"key": key, "size": size})
        trunc = root.find(f"{{{NS}}}IsTruncated")
        if trunc is not None and trunc.text == "true":
            continuation = root.find(f"{{{NS}}}NextContinuationToken").text
        else:
            break
    return objects


def parse_bam_filename(fname: str) -> dict:
    """
    Parse basecaller and modification info from BAM filename.
    Example: GM18501-ONT-hg38-R9-LSK110-guppy-sup-5mC.phased.bam
             GM19038-ONT-hg38-R10-LSK114-dorado081_sup_5mCG_5hmCG_v500.phased.bam
    """
    # Extract sample ID (everything before -ONT-)
    sample_match = re.match(r"^(GM\d+)-ONT-", fname)
    sample_id = sample_match.group(1) if sample_match else None

    # Pore chemistry
    pore = "R10" if "-R10-" in fname else "R9"

    # Modification calls present
    mods = []
    if "5hmCG" in fname or "5hmC" in fname:
        mods = ["5mC", "5hmC"]
    elif "5mCG" in fname:
        mods = ["5mC"]
    elif "5mC" in fname:
        mods = ["5mC"]

    # Basecaller
    if "guppy657" in fname:
        basecaller = "guppy657"
    elif "guppy" in fname:
        basecaller = "guppy"
    elif "dorado034" in fname:
        basecaller = "dorado034"
    elif "dorado050" in fname:
        basecaller = "dorado050"
    elif "dorado081" in fname:
        basecaller = "dorado081"
    elif "dorado082" in fname:
        basecaller = "dorado082"
    elif "dorado096" in fname:
        basecaller = "dorado096"
    else:
        basecaller = "unknown"

    # Has 5hmC?
    has_5hmc = "5hmC" in mods

    return {
        "sample_id": sample_id,
        "pore": pore,
        "basecaller": basecaller,
        "modifications": "+".join(mods),
        "has_5hmc": has_5hmc,
    }


import urllib.parse

def main():
    repo_root = Path(__file__).parents[2]
    out_path = repo_root / "config" / "bam_manifest.tsv"

    print("Listing BAMs from S3...", file=sys.stderr)
    objects = list_s3_prefix(BAM_PREFIX)
    bams = [o for o in objects if o["key"].endswith(".phased.bam")]
    print(f"Found {len(bams)} BAMs ({sum(o['size'] for o in bams)/1e12:.1f} TB)", file=sys.stderr)

    rows = []
    for bam in sorted(bams, key=lambda x: x["key"]):
        fname = Path(bam["key"]).name
        parsed = parse_bam_filename(fname)
        rows.append({
            "sample_id":     parsed["sample_id"],
            "bam_s3_url":    f"s3://1000g-ont/{bam['key']}",
            "bam_filename":  fname,
            "bam_size_gb":   round(bam["size"] / 1e9, 1),
            "pore":          parsed["pore"],
            "basecaller":    parsed["basecaller"],
            "modifications": parsed["modifications"],
            "has_5hmc":      str(parsed["has_5hmc"]).upper(),
        })

    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Written: {out_path}", file=sys.stderr)
    print(f"Samples: {len(rows)}", file=sys.stderr)

    # Summary stats
    r9 = sum(1 for r in rows if r["pore"] == "R9")
    r10 = sum(1 for r in rows if r["pore"] == "R10")
    with_5hmc = sum(1 for r in rows if r["has_5hmc"] == "TRUE")
    print(f"R9: {r9}  R10: {r10}  With 5hmC: {with_5hmc}  5mC-only: {len(rows)-with_5hmc}", file=sys.stderr)


if __name__ == "__main__":
    main()
