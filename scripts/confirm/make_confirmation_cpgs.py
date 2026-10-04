#!/usr/bin/env python3
"""make_confirmation_cpgs.py - regenerate config/confirmation_cpgs.txt.

Definition of 'prioritised hits' for the 1000-perm confirmation policy:
  best_cpg of the top-1000 prioritised variants (results/genome_452/)
  UNION the 90 SMR+HEIDI probes (eGFR 45 + T2D 23 + CAD 22).
Run from the repo root. The committed config/confirmation_cpgs.txt was
generated with exactly this script (1,084 CpGs, all 22 chromosomes).
"""
import collections
import pandas as pd

top = pd.read_csv("results/genome_452/top1000_prioritised_variants.tsv", sep="\t")
cpgs = set(top.best_cpg.dropna())
for f in ["eGFR_smr_heidi_hits.tsv", "t2d_smr_heidi_hits.tsv", "cvd_smr_heidi_hits.tsv"]:
    cpgs |= set(pd.read_csv(f"results/genome_452/coloc/{f}", sep="\t").probeID.dropna())
cpgs = sorted(cpgs)
with open("config/confirmation_cpgs.txt", "w") as o:
    o.write("\n".join(cpgs) + "\n")
cc = collections.Counter(c.split("_")[0] for c in cpgs)
print(f"wrote config/confirmation_cpgs.txt: {len(cpgs)} CpGs")
print(" ".join(f"{k}:{cc[k]}" for k in sorted(cc, key=lambda x: int(x.replace('chr', '')))))
