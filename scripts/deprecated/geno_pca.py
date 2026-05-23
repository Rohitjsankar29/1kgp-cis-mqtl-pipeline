#!/usr/bin/env python3
"""
DEPRECATED — replaced by plink2 --pca approx (Decision R16, 2026-03-14)

plink2 --set-all-var-ids '@:#:$r:$a' + --pca approx now handles genotype PCA
directly within 03_merge_vcfs_plink2.pbs. This script is no longer called by
any PBS job. Retained for reference only.

Original purpose: Compute genotype PCA from plink1 bed/bim/fam, writing
plink2-format eigenvec. Required at n=10 when plink2 --pca needed workarounds
for dot variant IDs. At n>=50, plink2 --pca approx works natively.
"""
import sys
import numpy as np
from pathlib import Path
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

def read_bed(prefix, n_samples, n_vars):
    bed_file = Path(prefix + '.bed')
    with open(bed_file, 'rb') as f:
        magic = f.read(3)
        assert magic == b'\x6c\x1b\x01', "Not a valid plink bed file (SNP-major mode required)"
        n_bytes_per_var = (n_samples + 3) // 4
        raw = np.frombuffer(f.read(), dtype=np.uint8)
    
    geno = np.full((n_vars, n_samples), np.nan)
    lookup = np.array([0, np.nan, 1, 2], dtype=float)  # 0=hom ref, 1=missing, 2=het, 3=hom alt
    
    for i in range(n_vars):
        block = raw[i * n_bytes_per_var:(i + 1) * n_bytes_per_var]
        bits = np.unpackbits(block, bitorder='little').reshape(-1, 2)[:n_samples]
        codes = bits[:, 0] + 2 * bits[:, 1]
        geno[i] = lookup[codes]
    
    return geno.T  # (n_samples, n_vars)

prefix = sys.argv[1]
n_pcs = int(sys.argv[2]) if len(sys.argv) > 2 else 10
out = sys.argv[3] if len(sys.argv) > 3 else prefix + '.eigenvec'

# Read fam
fam = []
with open(prefix + '.fam') as f:
    for line in f:
        parts = line.strip().split()
        fam.append((parts[0], parts[1]))  # FID, IID
n_samples = len(fam)
print(f"Samples: {n_samples}", flush=True)

# Count variants
with open(prefix + '.bim') as f:
    n_vars = sum(1 for _ in f)
print(f"Variants: {n_vars}", flush=True)

# Read genotype matrix
print("Reading bed file...", flush=True)
G = read_bed(prefix, n_samples, n_vars)

# Impute missing with column mean
print("Imputing missing genotypes...", flush=True)
col_means = np.nanmean(G, axis=0)
inds = np.where(np.isnan(G))
G[inds] = np.take(col_means, inds[1])

# Filter monomorphic variants
print("Filtering monomorphic variants...", flush=True)
var_mask = np.std(G, axis=0) > 0
G = G[:, var_mask]
print(f"Variants after monomorphic filter: {G.shape[1]}", flush=True)

# Standardise
print("Standardising...", flush=True)
G = StandardScaler().fit_transform(G)

# PCA
print(f"Running PCA (n_pcs={n_pcs})...", flush=True)
pca = PCA(n_components=n_pcs)
pcs = pca.fit_transform(G)
explained = pca.explained_variance_ratio_

print("Explained variance ratios:", explained, flush=True)

# Write eigenvec (plink2 format: #FID IID PC1 ... PCn)
with open(out, 'w') as f:
    header = '#FID\tIID\t' + '\t'.join(f'PC{i+1}' for i in range(n_pcs))
    f.write(header + '\n')
    for (fid, iid), row in zip(fam, pcs):
        vals = '\t'.join(f'{v:.6f}' for v in row)
        f.write(f'{fid}\t{iid}\t{vals}\n')

# Write eigenval
with open(out.replace('eigenvec', 'eigenval'), 'w') as f:
    for v in pca.explained_variance_:
        f.write(f'{v:.6f}\n')

print(f"Written: {out}", flush=True)
