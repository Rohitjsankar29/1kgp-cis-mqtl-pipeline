#!/usr/bin/env python3
"""Aggregate per-sample chr bedMethyl.gz files into tensorQTL phenotype BED.
Produces M-values and QC tables. Designed for pilot chromosomes first.
"""
import argparse, gzip, logging, sys
from pathlib import Path
import numpy as np, pandas as pd
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s', datefmt='%H:%M:%S')
log=logging.getLogger(__name__)

def read_one(path, min_cov):
    d={}
    with gzip.open(path,'rt') as f:
        for line in f:
            if line.startswith('#'): continue
            p=line.rstrip('\n').split('\t')
            if len(p)<11: continue
            cov=int(float(p[9])); pct=float(p[10])
            if cov>=min_cov: d[f'{p[0]}:{p[1]}']=pct/100.0
    return pd.Series(d, dtype=float)

def beta_to_m(x):
    x=np.clip(x,0.001,0.999)
    return np.log2(x/(1-x))

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--manifest', required=True)
    ap.add_argument('--bedmethyl-root', required=True)
    ap.add_argument('--chrom', default='chr22')
    ap.add_argument('--out-bed', required=True)
    ap.add_argument('--out-beta', default='')
    ap.add_argument('--out-pca', required=True)
    ap.add_argument('--out-qc-dir', required=True)
    ap.add_argument('--min-cov', type=int, default=5)
    ap.add_argument('--min-sample-frac', type=float, default=0.8)
    ap.add_argument('--n-pcs', type=int, default=3)
    args=ap.parse_args()
    samples=pd.read_csv(args.manifest, sep='\t')['sample_id'].tolist()
    root=Path(args.bedmethyl_root); qc=Path(args.out_qc_dir); qc.mkdir(parents=True, exist_ok=True)
    series={}; stats=[]
    for sid in samples:
        candidates=[root/args.chrom/sid/f'{sid}.{args.chrom}.cov{args.min_cov}.bedmethyl.gz', root/args.chrom/sid/f'{sid}.{args.chrom}.bedmethyl.gz']
        path=next((p for p in candidates if p.exists()), None)
        if path is None: sys.exit(f'Missing bedMethyl for {sid}: {candidates[0]}')
        log.info(f'Reading {sid}: {path}')
        s=read_one(path,args.min_cov); series[sid]=s
        stats.append({'sample_id':sid,'n_cpg_covered':len(s),'mean_beta':s.mean(),'median_beta':s.median()})
    pd.DataFrame(stats).to_csv(qc/'coverage_per_sample.tsv', sep='\t', index=False)
    beta=pd.DataFrame(series)
    valid=beta.notna().sum(axis=1)
    keep=valid>=int(np.ceil(args.min_sample_frac*len(samples)))
    pd.DataFrame({'cpg_id':beta.index,'n_valid':valid.values,'frac_valid':(valid/len(samples)).values,'passes_filter':keep.values}).to_csv(qc/'site_missingness.tsv', sep='\t', index=False)
    beta=beta.loc[keep]
    beta=beta.T.fillna(beta.mean(axis=1)).T
    if args.out_beta:
        write_bed(beta, args.out_beta, is_m=False)
    mval=beta.apply(beta_to_m)
    write_bed(mval, args.out_bed, is_m=True)
    mval.corr().to_csv(qc/'sample_correlation.tsv', sep='\t')
    from sklearn.decomposition import PCA
    n=min(args.n_pcs, mval.shape[1]-1, mval.shape[0])
    pcs=PCA(n_components=n).fit_transform(mval.T.values)
    pd.DataFrame(pcs.T, index=[f'MethPC{i+1}' for i in range(n)], columns=mval.columns).to_csv(args.out_pca, sep='\t')
    log.info(f'Wrote {args.out_bed}: {mval.shape[0]} CpGs x {mval.shape[1]} samples')

def write_bed(mat, out_path, is_m=True):
    chrom=[]; start=[]
    for x in mat.index:
        c,s=x.split(':'); chrom.append(c); start.append(int(s))
    bed=pd.DataFrame({'#chr':chrom,'start':start,'end':[s+1 for s in start],'phenotype_id':mat.index})
    bed=pd.concat([bed.reset_index(drop=True), mat.reset_index(drop=True)], axis=1).sort_values(['#chr','start'])
    with gzip.open(out_path,'wt') as f: bed.to_csv(f, sep='\t', index=False, float_format='%.6f')

if __name__=='__main__': main()
