#!/usr/bin/env python3
"""merge_confirmation.py - pool the per-chromosome 1000-perm confirmation
results, apply BH across the confirmation set, and report the confirmation
rate (overall + for the 90 SMR+HEIDI probes).

Outputs:
  {out_dir}/confirmation_1000perm.tsv       per-CpG confirmation table (+ qval)
  {out_dir}/confirmation_summary.txt        headline confirmation numbers
"""
import argparse, glob
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--in-glob",
                default="/g/data/cy94/rs4477/downstream/genome_452/chr*/confirm/*.confirm1000.txt.gz")
ap.add_argument("--smr-dir", default=None,
                help="dir with {eGFR,t2d,cvd}_smr_heidi_hits.tsv to cross-mark")
ap.add_argument("--out-dir", default="/g/data/cy94/rs4477/downstream/genome_452")
a = ap.parse_args()

files = sorted(glob.glob(a.in_glob))
if not files:
    raise SystemExit(f"no files match {a.in_glob}")
d = pd.concat([pd.read_csv(f, sep="\t") for f in files], ignore_index=True)
pcol = "pval_beta" if "pval_beta" in d.columns else \
       ("pval_perm" if "pval_perm" in d.columns else "pval_nominal")
print(f"pooled {len(d):,} CpGs from {len(files)} files; p-value column: {pcol}")

p = d[pcol].to_numpy(float)
order = np.argsort(p)
n = len(p)
q = np.empty(n)
q[order] = np.minimum.accumulate((p[order] * n / np.arange(1, n + 1))[::-1])[::-1]
d["qval_confirm"] = np.clip(q, 0, 1)
d["confirmed_q05"] = d["qval_confirm"] < 0.05
d.sort_values(pcol).to_csv(f"{a.out_dir}/confirmation_1000perm.tsv", sep="\t", index=False)

idcol = "phenotype_id" if "phenotype_id" in d.columns else d.columns[0]
lines = ["1000-PERMUTATION CONFIRMATION - prioritised hits (top1000 + SMR probes)",
         "=" * 68,
         f"CpGs tested: {len(d):,}",
         f"confirmed (BH q<0.05 within confirmation set): {int(d.confirmed_q05.sum()):,} "
         f"({100 * d.confirmed_q05.mean():.1f}%)"]
if a.smr_dir:
    smr = set()
    for f in ["eGFR_smr_heidi_hits.tsv", "t2d_smr_heidi_hits.tsv", "cvd_smr_heidi_hits.tsv"]:
        try:
            smr |= set(pd.read_csv(f"{a.smr_dir}/{f}", sep="\t").probeID)
        except FileNotFoundError:
            pass
    m = d[idcol].isin(smr)
    lines.append(f"SMR+HEIDI probes present: {int(m.sum())} / {len(smr)}; "
                 f"confirmed q<0.05: {int(d.loc[m, 'confirmed_q05'].sum())}")
with open(f"{a.out_dir}/confirmation_summary.txt", "w") as o:
    o.write("\n".join(lines) + "\n")
print("\n".join(lines))
