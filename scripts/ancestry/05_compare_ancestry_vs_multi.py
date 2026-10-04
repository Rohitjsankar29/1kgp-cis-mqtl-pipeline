#!/usr/bin/env python3
"""05_compare_ancestry_vs_multi.py — the payoff script.

For each SMR+HEIDI-significant locus, join the primary (multi-ancestry LD)
result with the ancestry-restricted re-run and classify the locus into one of
three categories:

    ROBUST      — significant under both  (mediation call unchanged)
    LOST        — significant under multi-ancestry, not under the ancestry-restricted run
                  (candidate false positive driven by LD mismatch)
    GAINED      — the opposite  (candidate false negative under the primary run;
                  rare when the restricted panel is smaller than the primary)

Also joins the coloc.susie PP-vector where present, so a locus can be
"SMR-ROBUST but coloc-DISCORDANT" — flagged separately.

Writes:
    ancestry_sensitivity_<POP>_<TRAIT>.tsv   (per-locus classified table)
    ancestry_sensitivity_summary.tsv         (counts across pop x trait)

Usage:
    python 05_compare_ancestry_vs_multi.py \
        --pop EUR \
        --primary-dir  /g/data/cy94/rs4477/downstream/genome_452/smr \
        --ancestry-dir /g/data/cy94/rs4477/downstream/genome_452/ancestry/smr/EUR \
        --coloc-primary  /g/data/cy94/rs4477/downstream/genome_452/coloc \
        --coloc-ancestry /g/data/cy94/rs4477/downstream/genome_452/ancestry/coloc/EUR \
        --out-dir /g/data/cy94/rs4477/downstream/genome_452/ancestry/sensitivity
"""
import argparse, glob, os
import pandas as pd

TRAITS = ["egfr", "cad", "t2d"]
SMR_P_THRESH = 5e-8 / 90   # Bonferroni-adjusted for the 90 candidate probes
HEIDI_P_THRESH = 0.05


def load_smr(path_pattern):
    files = glob.glob(path_pattern)
    if not files:
        return pd.DataFrame()
    frames = [pd.read_csv(f, sep="\t") for f in files]
    return pd.concat(frames, ignore_index=True)


def classify(primary_sig, ancestry_sig):
    """Three-way classification of a locus."""
    if primary_sig and ancestry_sig:      return "ROBUST"
    if primary_sig and not ancestry_sig:  return "LOST"
    if not primary_sig and ancestry_sig:  return "GAINED"
    return "NEITHER"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pop", required=True)
    ap.add_argument("--primary-dir",  required=True)
    ap.add_argument("--ancestry-dir", required=True)
    ap.add_argument("--coloc-primary",  default=None)
    ap.add_argument("--coloc-ancestry", default=None)
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()

    os.makedirs(a.out_dir, exist_ok=True)
    summary_rows = []

    for trait in TRAITS:
        # SMR outputs — filename convention: smr_<trait>.smr
        prim = load_smr(os.path.join(a.primary_dir,  f"smr_{trait}*.smr"))
        anc  = load_smr(os.path.join(a.ancestry_dir, f"smr_{trait}*.smr"))
        if prim.empty or anc.empty:
            print(f"{trait}: missing SMR outputs (primary or ancestry) — skipping")
            continue

        # SMR-significant AND HEIDI-passing
        def sig(df):
            return (df["p_SMR"] < SMR_P_THRESH) & (df["p_HEIDI"].fillna(1.0) > HEIDI_P_THRESH)
        prim["sig_primary"]   = sig(prim)
        anc["sig_ancestry"]   = sig(anc)

        m = prim.merge(anc, on="probeID", suffixes=("_primary", "_ancestry"), how="outer")
        m[["sig_primary", "sig_ancestry"]] = m[["sig_primary", "sig_ancestry"]].fillna(False)
        m["class"] = m.apply(lambda r: classify(r["sig_primary"], r["sig_ancestry"]), axis=1)

        # optional coloc merge
        if a.coloc_primary and a.coloc_ancestry:
            cp = load_smr(os.path.join(a.coloc_primary, trait, "coloc_susie_*.tsv"))
            ca = load_smr(os.path.join(a.coloc_ancestry, trait, "coloc_susie_*.tsv"))
            if not cp.empty and "PP.H4.abf" in cp.columns:
                cp = cp.rename(columns={"PP.H4.abf": "PP4_primary"})[["locus", "PP4_primary"]]
                ca = ca.rename(columns={"PP.H4.abf": "PP4_ancestry"})[["locus", "PP4_ancestry"]]
                m = m.merge(cp, left_on="probeID", right_on="locus", how="left")
                m = m.merge(ca, left_on="probeID", right_on="locus", how="left")
                m["coloc_concordant"] = ((m["PP4_primary"] > 0.8) == (m["PP4_ancestry"] > 0.8))

        outf = os.path.join(a.out_dir, f"ancestry_sensitivity_{a.pop}_{trait}.tsv")
        m.to_csv(outf, sep="\t", index=False)

        counts = m["class"].value_counts().to_dict()
        summary_rows.append({
            "pop": a.pop, "trait": trait,
            "n_probes": len(m),
            "n_ROBUST": counts.get("ROBUST", 0),
            "n_LOST":   counts.get("LOST", 0),
            "n_GAINED": counts.get("GAINED", 0),
        })
        print(f"{a.pop} · {trait}: {counts}")

    summary = pd.DataFrame(summary_rows)
    summary_out = os.path.join(a.out_dir, f"ancestry_sensitivity_summary_{a.pop}.tsv")
    summary.to_csv(summary_out, sep="\t", index=False)
    print(f"\nsummary → {summary_out}")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
