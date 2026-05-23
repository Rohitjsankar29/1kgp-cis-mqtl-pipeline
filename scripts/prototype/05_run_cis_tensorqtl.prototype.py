#!/usr/bin/env python3
"""Run TensorQTL cis mapping with optional phenotype chunking."""
import argparse, logging, os, sys
import pandas as pd
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s', datefmt='%H:%M:%S')
log=logging.getLogger(__name__)

def load_pheno(path,n_chunks,chunk_idx):
    df=pd.read_csv(path, sep='\t')
    df=df.sort_values(['#chr','start']).set_index('phenotype_id')
    pos=df[['#chr','end']].rename(columns={'#chr':'chr','end':'pos'})
    ph=df.drop(columns=['#chr','start','end']).astype(float)
    if n_chunks>1:
        total=len(ph); size=(total+n_chunks-1)//n_chunks; a=chunk_idx*size; b=min(a+size,total)
        ph=ph.iloc[a:b]; pos=pos.iloc[a:b]
        log.info(f'Phenotype chunk {chunk_idx}/{n_chunks}: {a}-{b}')
    return ph,pos

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--plink-prefix', required=True)
    ap.add_argument('--phenotype-bed', required=True)
    ap.add_argument('--covariates', required=True)
    ap.add_argument('--out-prefix', required=True)
    ap.add_argument('--cis-window', type=int, default=1000000)
    ap.add_argument('--n-permutations', type=int, default=1000)
    ap.add_argument('--n-chunks', type=int, default=1)
    ap.add_argument('--chunk-idx', type=int, default=0)
    ap.add_argument('--nominal', action='store_true')
    args=ap.parse_args()
    import torch, tensorqtl
    from tensorqtl import cis, genotypeio
    log.info(f'CUDA: {torch.cuda.is_available()} tensorQTL: {tensorqtl.__version__}')
    pr=genotypeio.PlinkReader(args.plink_prefix)
    bim=pr.bim.copy(); bim['vid']=bim['snp']
    geno=pr.load_genotypes(); geno.index=bim['vid'].values
    var=bim.set_index('vid')[['chrom','pos']].rename(columns={'chrom':'chr'})
    ph,pos=load_pheno(args.phenotype_bed,args.n_chunks,args.chunk_idx)
    cov_raw=pd.read_csv(args.covariates, sep='\t', index_col=0)
    common=sorted(set(pr.fam['iid']) & set(ph.columns) & set(cov_raw.columns))
    if len(common)<5: sys.exit(f'Too few common samples: {len(common)}')
    geno=geno[common]; ph=ph[common]; cov=cov_raw[common].T
    os.makedirs(os.path.dirname(args.out_prefix) or '.', exist_ok=True)
    if args.nominal:
        cis.map_nominal(geno,var,ph,pos,args.out_prefix,covariates_df=cov,window=args.cis_window)
    perm=cis.map_cis(geno,var,ph,pos,covariates_df=cov,nperm=args.n_permutations,window=args.cis_window)
    perm.to_csv(args.out_prefix+'.permutation.txt.gz', sep='\t', compression='gzip')
    log.info('Done')
if __name__=='__main__': main()
