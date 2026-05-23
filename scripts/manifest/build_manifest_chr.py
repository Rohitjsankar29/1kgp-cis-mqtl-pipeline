"""
build_manifest_chr.py
─────────────────────
Build config/bam_manifest.tsv from the 1000G-ONT S3 bucket, then run
modkit pileup for a single user-specified chromosome across all samples.

Usage
-----
    python build_manifest_chr.py <chromosome>

Examples
--------
    python build_manifest_chr.py chr22
    python build_manifest_chr.py chr1
    python build_manifest_chr.py chrX

Requirements (must be on PATH)
-------------------------------
    modkit   tabix   bgzip   samtools

Output layout
-------------
    config/bam_manifest.tsv
    output/per_chr/<chrom>/<sample_id>/
        <sample_id>.<chrom>.bedmethyl.gz        (all CpGs)
        <sample_id>.<chrom>.bedmethyl.gz.tbi
        <sample_id>.<chrom>.bedmethyl.cov10.gz  (>=10x coverage)
        <sample_id>.<chrom>.bedmethyl.cov10.gz.tbi
        <sample_id>.<chrom>.log
"""

# ── standard library ────────────────────────────────────────────────────────
import csv
import re
import subprocess
import sys
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

# ════════════════════════════════════════════════════════════════════════════
#  USER CONFIGURATION  ── edit these two lines before running
# ════════════════════════════════════════════════════════════════════════════
REF    = "/path/to/hg38.fa"       # <-- UPDATE THIS: absolute path to hg38 FASTA
OUTDIR = Path("output/per_chr")   # <-- UPDATE THIS: where results are written
# ════════════════════════════════════════════════════════════════════════════

# S3 / API constants
S3_BASE    = "https://s3.amazonaws.com/1000g-ont"
BAM_PREFIX = "PROCESSED_DATA/ALIGNED_TO_HG38/MINIMAP2_ALIGNED_BAMS/"
NS         = "http://s3.amazonaws.com/doc/2006-03-01/"

# Runtime defaults
THREADS    = 16
COV_CUTOFF = 10   # minimum read depth to keep a CpG site


# ── helpers ──────────────────────────────────────────────────────────────────

def check_tools() -> None:
    """Exit early if any required tool is missing from PATH."""
    for tool in ("modkit", "tabix", "bgzip", "samtools"):
        result = subprocess.run(
            ["which", tool], capture_output=True, text=True
        )
        if result.returncode != 0:
            print(f"ERROR: '{tool}' not found on PATH. Please install it first.",
                  file=sys.stderr)
            sys.exit(1)
    print("All required tools found.", file=sys.stderr)


def list_s3_prefix(prefix: str) -> list[dict]:
    """
    Return every object under *prefix* in the S3 bucket.
    Handles pagination automatically (S3 returns max 1 000 keys per request).
    """
    objects      = []
    continuation = None

    while True:
        url = (
            f"{S3_BASE}?list-type=2"
            f"&prefix={urllib.parse.quote(prefix)}"
            f"&max-keys=1000"
        )
        if continuation:
            url += f"&continuation-token={urllib.parse.quote(continuation)}"

        with urllib.request.urlopen(url) as response:
            root = ET.parse(response).getroot()

        for item in root.findall(f"{{{NS}}}Contents"):
            objects.append({
                "key":  item.find(f"{{{NS}}}Key").text,
                "size": int(item.find(f"{{{NS}}}Size").text),
            })

        truncated = root.find(f"{{{NS}}}IsTruncated")
        if truncated is not None and truncated.text == "true":
            continuation = root.find(f"{{{NS}}}NextContinuationToken").text
        else:
            break   # no more pages

    return objects


def parse_bam_filename(fname: str) -> dict:
    """
    Extract sample_id, pore chemistry, basecaller, and modification
    model directly from the BAM filename.

    Supported filename patterns
    ───────────────────────────
    GM18501-ONT-hg38-R9-LSK110-guppy-sup-5mC.phased.bam
    GM19038-ONT-hg38-R10-LSK114-dorado081_sup_5mCG_5hmCG_v500.phased.bam
    """
    # Sample ID: everything before the first -ONT-
    sample_match = re.match(r"^(GM\d+)-ONT-", fname)
    sample_id    = sample_match.group(1) if sample_match else "unknown"

    # Pore chemistry
    pore = "R10" if "-R10-" in fname else "R9"

    # Modification model (order matters: check 5hmC before plain 5mC)
    if "5hmCG" in fname or "5hmC" in fname:
        mods = ["5mC", "5hmC"]
    elif "5mCG" in fname or "5mC" in fname:
        mods = ["5mC"]
    else:
        mods = []

    # Basecaller — regex captures any dorado version number automatically
    dorado_match = re.search(r"dorado(\d+)", fname)
    if "guppy657" in fname:
        basecaller = "guppy657"
    elif "guppy" in fname:
        basecaller = "guppy"
    elif dorado_match:
        basecaller = f"dorado{dorado_match.group(1)}"
    else:
        basecaller = "unknown"

    return {
        "sample_id":     sample_id,
        "pore":          pore,
        "basecaller":    basecaller,
        "modifications": "+".join(mods),
        "has_5hmc":      "5hmC" in mods,
    }


# ── core analysis ────────────────────────────────────────────────────────────

def run_modkit_for_chr(
    sample_id:   str,
    bam_s3_url:  str,
    chrom:       str,
    mods:        str,
) -> None:
    """
    Stream a single chromosome from a remote S3 BAM with samtools, pipe
    directly into modkit pileup, then filter to >=COV_CUTOFF coverage.

    The '-' argument to modkit tells it to read from stdin, so no
    intermediate BAM file is written to disk.
    """
    out_dir = OUTDIR / chrom / sample_id
    out_dir.mkdir(parents=True, exist_ok=True)

    outfile  = out_dir / f"{sample_id}.{chrom}.bedmethyl.gz"
    filtered = out_dir / f"{sample_id}.{chrom}.bedmethyl.cov{COV_CUTOFF}.gz"
    logfile  = out_dir / f"{sample_id}.{chrom}.log"

    # Choose the right modification flag for modkit
    mod_flag = "5mCG_5hmCG" if "5hmC" in mods else "5mC"

    print(f"  [{sample_id}] {chrom} — running modkit pileup …", file=sys.stderr)

    # Step 1 ── stream chromosome from S3 and run modkit
    step1 = (
        f"samtools view -b -@ {THREADS} {bam_s3_url} {chrom} | "
        f"modkit pileup "
        f"  --cpg "
        f"  --ref {REF} "
        f"  --modified-bases {mod_flag} "
        f"  --combine-strands "
        f"  --threads {THREADS} "
        f"  --bgzf "
        f"  --log {logfile} "
        f"  - {outfile}"
    )
    result = subprocess.run(step1, shell=True)
    if result.returncode != 0:
        print(f"  ERROR: modkit failed for {sample_id} on {chrom} — skipping.",
              file=sys.stderr)
        return

    # Step 2 ── index the raw output
    subprocess.run(f"tabix -p bed {outfile}", shell=True, check=True)

    # Step 3 ── coverage filter (bgzip -dc is safer than zcat for bgzf files)
    subprocess.run(
        f"bgzip -dc {outfile} | awk '$10 >= {COV_CUTOFF}' | bgzip > {filtered}",
        shell=True, check=True,
    )

    # Step 4 ── index the filtered output
    subprocess.run(f"tabix -p bed {filtered}", shell=True, check=True)

    # Step 5 ── quick line counts to log
    total = int(subprocess.check_output(
        f"bgzip -dc {outfile} | wc -l", shell=True
    ).strip())
    kept = int(subprocess.check_output(
        f"bgzip -dc {filtered} | wc -l", shell=True
    ).strip())
    print(
        f"  [{sample_id}] {chrom} done — "
        f"total CpGs: {total:,}  kept (>={COV_CUTOFF}x): {kept:,}",
        file=sys.stderr,
    )


# ── manifest builder ─────────────────────────────────────────────────────────

def build_manifest(bams: list[dict]) -> list[dict]:
    """Parse all BAM entries and return manifest rows."""
    rows = []
    for bam in sorted(bams, key=lambda x: x["key"]):
        fname  = Path(bam["key"]).name
        parsed = parse_bam_filename(fname)
        rows.append({
            "sample_id":     parsed["sample_id"],
            "bam_s3_url":    f"s3://1000g-ont/{bam['key']}",
            "bam_filename":  fname,
            "bam_size_gb":   round(bam["size"] / 1e9, 1),
            "pore":          parsed["pore"],
            "basecaller":    parsed["basecaller"],
            "modifications": parsed["modifications"],
            "has_5hmc":      "TRUE" if parsed["has_5hmc"] else "FALSE",
        })
    return rows


def write_manifest(rows: list[dict], out_path: Path) -> None:
    """Write manifest rows to a tab-separated file."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Manifest written → {out_path}", file=sys.stderr)


def print_summary(rows: list[dict]) -> None:
    """Print a brief summary of the manifest to stderr."""
    r9        = sum(1 for r in rows if r["pore"] == "R9")
    r10       = sum(1 for r in rows if r["pore"] == "R10")
    with_5hmc = sum(1 for r in rows if r["has_5hmc"] == "TRUE")
    print(
        f"Samples: {len(rows)}  |  R9: {r9}  R10: {r10}  "
        f"5mC+5hmC: {with_5hmc}  5mC-only: {len(rows) - with_5hmc}",
        file=sys.stderr,
    )


# ── entry point ──────────────────────────────────────────────────────────────

def main() -> None:
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    chrom = sys.argv[1]
    print(f"Target chromosome : {chrom}", file=sys.stderr)
    print(f"Reference         : {REF}",   file=sys.stderr)
    print(f"Output directory  : {OUTDIR}", file=sys.stderr)

    # Pre-flight: confirm all tools are available
    check_tools()

    # ── 1. Build manifest ────────────────────────────────────────────────────
    print("\nListing BAMs from S3 …", file=sys.stderr)
    objects = list_s3_prefix(BAM_PREFIX)
    bams    = [o for o in objects if o["key"].endswith(".phased.bam")]
    total_tb = sum(o["size"] for o in bams) / 1e12
    print(f"Found {len(bams)} BAMs  ({total_tb:.1f} TB total)", file=sys.stderr)

    rows = build_manifest(bams)

    repo_root    = Path(__file__).parents[2]
    manifest_out = repo_root / "config" / "bam_manifest.tsv"
    write_manifest(rows, manifest_out)
    print_summary(rows)

    # ── 2. Run modkit for the requested chromosome ───────────────────────────
    print(f"\nRunning modkit pileup for {chrom} on {len(rows)} samples …\n",
          file=sys.stderr)

    for row in rows:
        run_modkit_for_chr(
            sample_id  = row["sample_id"],
            bam_s3_url = row["bam_s3_url"],
            chrom      = chrom,
            mods       = row["modifications"],
        )

    print(f"\nAll samples complete for {chrom}.", file=sys.stderr)


if __name__ == "__main__":
    main()
