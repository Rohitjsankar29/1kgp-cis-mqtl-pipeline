#!/usr/bin/env python3
"""Build a simple CpG-centric variant prioritisation feature matrix."""
import argparse, numpy as np, pandas as pd

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--qtl', required=True, help='TensorQTL nominal/permutation-like TSV')
    ap.add_argument('--out', required=True)
    args=ap.parse_args()
    df=pd.read_csv(args.qtl, sep='\t')
    if 'phenotype_id' not in df.columns and df.index.name: df=df.reset_index().rename(columns={df.index.name:'phenotype_id'})
    if 'pval_nominal' in df.columns: p='pval_nominal'
    elif 'pval_beta' in df.columns: p='pval_beta'
    else: raise SystemExit('Need pval_nominal or pval_beta')
    if 'slope' not in df.columns: df['slope']=df.get('beta_shape1',0)
    df['neg_log10_p']=-np.log10(df[p].clip(1e-300))
    df['abs_effect']=df['slope'].abs()
    if 'variant_id' in df.columns:
        df['variant_pos']=df.variant_id.astype(str).str.split(':').str[1].astype(float)
    else: df['variant_pos']=np.nan
    df['cpg_pos']=df.phenotype_id.astype(str).str.split(':').str[1].astype(float)
    df['distance_to_cpg']=(df.variant_pos-df.cpg_pos).abs()
    # placeholders for later annotation/finemapping/SV integration
    df['pip']=df.get('pip',0)
    df['is_sv']=df.get('is_sv',0)
    df['annotation_score']=df.get('annotation_score',0)
    feats=['neg_log10_p','abs_effect','pip','is_sv','annotation_score']
    for c in feats+['distance_to_cpg']:
        x=df[c].astype(float); df[c+'_z']=(x-x.mean())/(x.std(ddof=0)+1e-9)
    df['distance_score_z']=-df['distance_to_cpg_z']
    df['priority_score']=df['neg_log10_p_z']+df['abs_effect_z']+df['pip_z']+df['is_sv_z']+df['annotation_score_z']+df['distance_score_z']
    df.sort_values(['phenotype_id','priority_score'], ascending=[True,False]).to_csv(args.out, sep='\t', index=False)
if __name__=='__main__': main()
