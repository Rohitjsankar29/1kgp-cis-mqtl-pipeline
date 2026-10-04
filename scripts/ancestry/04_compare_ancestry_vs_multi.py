#!/usr/bin/env python3
"""04_compare_ancestry_vs_multi.py -- the payoff script.

For each SMR probe (CpG) tested against each of the three cardiometabolic
GWAS, join the PRIMARY (multi-ancestry LD, n=452) SMR output with the
ANCESTRY-RESTRICTED (e.g. EUR n=71) SMR output and classify each probe:

    ROBUST      -- SMR-significant AND HEIDI-passing under both LD references
    LOST        -- significant under multi-ancestry only  (candidate FP from LD mismatch)
    GAINED      -- significant under the ancestry LD only (rare with smaller panels)
    NEITHER     -- non-significant under both

Writes per-trait per-locus tables and a per-superpop summary. This is the
deliverable that lets the presentation say "we tested this and here's what
we found" instead of leaving the ancestry limitation as a footnote.

Usage:
    python3 04_compare_ancestry_vs_multi.py \\
        --pop EUR \\
        --primary-dir  /scratch/cy94/rs4477/coloc \\
        --ancestry-dir /g/data/cy94/rs4477/downstream/genome_452/ancestry/smr/EUR \\
        --out-dir      /g/data/cy94/rs4477/downstream/genome_452/ancestry/sensitivity
"""
import argparse, os, sys
import pandas as pd

TRAITS = ["eGFR", "cvd", "t2d"]
HEIDI_P_THRESH = 0.05     # HEIDI p > 0.05 = pass (consistent with mediation)


def load_smr_all(directory, trait):
    """Load the concatenated 'smr_<trait>_all.smr' file."""
    path = os.path.join(directory, f"smr_{trait}_all.smr")
    if not os.path.isfile(path):
        # try per-chromosome fallback
        parts = []
        for c in range(1, 23):
            p = os.path.join(directory, f"smr_{trait}_chr{c}.smr")
            if os.path.isfile(p):
                parts.append(pd.read_csv(p, sep="\t"))
        if not parts:
            print(f"WARNING: no SMR output under {directory} for {trait}", file=sys.stderr)
            return None
        return pd.concat(parts, ignore_index=True)
    return pd.read_csv(path, sep="\t")


def sig_mask(df, tot):
    """SMR-significant (Bonferroni over the tested set) AND HEIDI-passing."""
    p_smr = df.get("p_SMR", df.get("pSMR"))
    p_heidi = df.get("p_HEIDI", df.get("pHEIDI"))
    if p_smr is None or p_heidi is None:
        raise SystemExit(f"missing p_SMR or p_HEIDI columns; have {list(df.columns)}")
    return (p_smr.astype(float) < 0.05 / tot) & (p_heidi.fillna(1.0).astype(float) > HEIDI_P_THRESH)


def classify(prim, anc):
    if prim and anc:      return "ROBUST"
    if prim and not anc:  return "LOST"
    if not prim and anc:  return "GAINED"
    return "NEITHER"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pop", required=True)
    ap.add_argument("--primary-dir",  required=True,
                    help="e.g. /scratch/cy94/rs4477/coloc  (has smr_<trait>_all.smr)")
    ap.add_argument("--ancestry-dir", required=True,
                    help="e.g. .../ancestry/smr/EUR")
    ap.add_argument("--out-dir",  required=True)
    a = ap.parse_args()

    os.makedirs(a.out_dir, exist_ok=True)
    summary_rows = []

    for trait in TRAITS:
        prim = load_smr_all(a.primary_dir,  trait)
        anc  = load_smr_all(a.ancestry_dir, trait)
        if prim is None or anc is None:
            print(f"{trait}: skipping (missing input)")
            continue

        # Bonferroni denominator = number tested in the primary set
        tot = len(prim)
        prim["sig_primary"]  = sig_mask(prim, tot).values
        anc["sig_ancestry"]  = sig_mask(anc,  tot).values

        key = "probeID" if "probeID" in prim.columns else prim.columns[0]
        m = prim.merge(anc, on=key, suffixes=("_primary", "_ancestry"), how="outer")
        m["sig_primary"]  = m["sig_primary"].fillna(False)
        m["sig_ancestry"] = m["sig_ancestry"].fillna(False)
        m["class"] = m.apply(lambda r: classify(r["sig_primary"], r["sig_ancestry"]), axis=1)

        outf = os.path.join(a.out_dir, f"ancestry_sensitivity_{a.pop}_{trait}.tsv")
        m.to_csv(outf, sep="\t", index=False)

        counts = m["class"].value_counts().to_dict()
        # only report probes flagged by at least one LD reference
        counts_flagged = {k: v for k, v in counts.items() if k != "NEITHER"}
        summary_rows.append({
            "pop": a.pop, "trait": trait, "n_probes_tested": tot,
            "n_ROBUST": counts.get("ROBUST", 0),
            "n_LOST":   counts.get("LOST", 0),
            "n_GAINED": counts.get("GAINED", 0),
        })
        print(f"{a.pop} · {trait}:  ROBUST={counts.get('ROBUST',0)}  "
              f"LOST={counts.get('LOST',0)}  GAINED={counts.get('GAINED',0)}  "
              f"(of {tot} tested)")

    summary = pd.DataFrame(summary_rows)
    summary_out = os.path.join(a.out_dir, f"ancestry_sensitivity_summary_{a.pop}.tsv")
    summary.to_csv(summary_out, sep="\t", index=False)
    print(f"\nsummary -> {summary_out}\n{summary.to_string(index=False)}")


if __name__ == "__main__":
    main()
