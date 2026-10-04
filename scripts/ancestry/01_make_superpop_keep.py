#!/usr/bin/env python3
"""01_make_superpop_keep.py -- build PLINK --keep files per superpopulation.

Reads config/metadata_with_pop.tsv (produced by 00_prepare_superpop_labels.sh)
and writes keep_<POP>.txt lists in the requested output directory. Each list
uses PLINK's FID<TAB>IID format; for 1000G samples the FID equals the IID.

Usage:
    python3 01_make_superpop_keep.py \\
        --metadata /scratch/cy94/rs4477/1kgp-cis-mqtl-pipeline/config/metadata_with_pop.tsv \\
        --out-dir  /g/data/cy94/rs4477/downstream/genome_452/ancestry/keep
"""
import argparse, os, sys
import pandas as pd

SUPERPOPS = ["EUR", "AFR", "SAS", "EAS", "AMR"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--out-dir",  required=True)
    ap.add_argument("--iid-col",  default="sample",
                    help="column name for sample/IID (default: sample)")
    ap.add_argument("--pop-col",  default="superpopulation",
                    help="column name for superpopulation code (default: superpopulation)")
    a = ap.parse_args()

    os.makedirs(a.out_dir, exist_ok=True)
    m = pd.read_csv(a.metadata, sep="\t")
    for c in (a.iid_col, a.pop_col):
        if c not in m.columns:
            sys.exit(f"metadata missing column {c!r}; have {list(m.columns)}")

    for pop in SUPERPOPS:
        sub = m[m[a.pop_col] == pop]
        path = os.path.join(a.out_dir, f"keep_{pop}.txt")
        # PLINK --keep expects FID<TAB>IID; 1000G FID == IID
        sub[[a.iid_col]].assign(FID="0").reindex(columns=["FID", a.iid_col]).to_csv(path, sep="\t", header=False, index=False)
        print(f"{pop}\t{len(sub)}\t{path}")


if __name__ == "__main__":
    main()
