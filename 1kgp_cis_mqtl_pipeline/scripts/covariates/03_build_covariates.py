#!/usr/bin/env python3
"""Build tensorQTL covariates: rows=covariates, columns=samples."""
import argparse, logging, sys
import pandas as pd
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s', datefmt='%H:%M:%S')
log=logging.getLogger(__name__)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--manifest', required=True)
    ap.add_argument('--meth-pca', required=True)
    ap.add_argument('--geno-pca', default='')
    ap.add_argument('--n-meth-pcs', type=int, default=3)
    ap.add_argument('--n-geno-pcs', type=int, default=3)
    ap.add_argument('--out', required=True)
    args=ap.parse_args()
    m=pd.read_csv(args.manifest, sep='\t')
    ids=m.sample_id.tolist(); rows={}
    if 'sex' in m.columns:
        smap={'male':1,'m':1,'xy':1,'female':0,'f':0,'xx':0,'NA':0}
        rows['sex_male']={r.sample_id:smap.get(str(r.sex).lower(),0) for _,r in m.iterrows()}
    if 'basecaller_model' in m.columns:
        rows['basecaller_model']={r.sample_id:r.basecaller_model for _,r in m.iterrows()}
    if args.geno_pca:
        g=pd.read_csv(args.geno_pca, sep=r'\s+')
        iid='IID' if 'IID' in g.columns else g.columns[1]
        g=g.set_index(iid)
        for i in range(1,args.n_geno_pcs+1):
            pc=f'PC{i}'
            if pc in g.columns: rows[f'GenoPC{i}']={sid:g.loc[sid,pc] for sid in ids if sid in g.index}
    mp=pd.read_csv(args.meth_pca, sep='\t', index_col=0)
    for i in range(1,args.n_meth_pcs+1):
        pc=f'MethPC{i}'
        if pc in mp.index: rows[pc]={sid:mp.loc[pc,sid] for sid in ids if sid in mp.columns}
    cov=pd.DataFrame(rows).T.reindex(columns=ids).fillna(0)
    cov.to_csv(args.out, sep='\t', float_format='%.6f')
    log.info(f'Wrote {args.out}: {cov.shape[0]} covariates x {cov.shape[1]} samples')
if __name__=='__main__': main()
