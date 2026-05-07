#!/usr/bin/env python3
"""Small QC for TensorQTL permutation results."""
import argparse, os, numpy as np, pandas as pd
from scipy.stats import chi2
from statsmodels.stats.multitest import multipletests

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--permutation', required=True); ap.add_argument('--out-dir', required=True); ap.add_argument('--label', default='pilot'); ap.add_argument('--fdr', type=float, default=0.05); args=ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    df=pd.read_csv(args.permutation, sep='\t', index_col=0)
    p=df['pval_beta'].dropna().clip(1e-300,1-1e-16)
    lam=np.median(chi2.ppf(1-p,1))/chi2.ppf(0.5,1) if len(p)>10 else np.nan
    q=multipletests(p, method='fdr_bh')[1]
    sig=df.loc[p.index].copy(); sig['qval']=q; sig=sig[sig.qval<args.fdr]
    sig.to_csv(f'{args.out_dir}/{args.label}.significant.tsv', sep='\t')
    open(f'{args.out_dir}/{args.label}.lambda_gc.txt','w').write(f'label\tlambda_gc\tn_tests\n{args.label}\t{lam:.4f}\t{len(p)}\n')
    print(f'lambda={lam:.4f}; significant={len(sig)}')
if __name__=='__main__': main()
