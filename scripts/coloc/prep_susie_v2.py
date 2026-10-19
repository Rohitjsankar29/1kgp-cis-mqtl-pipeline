import pandas as pd, numpy as np, os, pyarrow.parquet as pq
from pandas_plink import read_plink1_bin
G="/g/data/cy94/rs4477/downstream/genome_452"
OUT="/scratch/cy94/rs4477/coloc/susie_loci"; os.makedirs(OUT,exist_ok=True)
hits=pd.read_csv("/scratch/cy94/rs4477/coloc/eGFR_smr_heidi_hits.tsv",sep="\t")
gw=pd.read_csv("/scratch/cy94/rs4477/coloc/gwas/eGFR_hg38.ma.full",sep="\t")
gwi=gw.drop_duplicates('SNP').set_index('SNP')
n_ok=0
plink_cache={}
for _,hit in hits.iterrows():
    cpg=hit.probeID; chrom=cpg.split('_')[0]
    pth=f"{G}/{chrom}/nominal/{chrom}.sig.cis_qtl_pairs.{chrom}.parquet"
    if chrom in ['chr1','chr2']: pth=f"{G}/{chrom}/nominal/{chrom}.sig.merged.parquet"
    pf=pq.ParquetFile(pth); fr=[]
    for i in range(pf.num_row_groups):
        t=pf.read_row_group(i,columns=['phenotype_id','variant_id','af','slope','slope_se']).to_pandas()
        fr.append(t[t.phenotype_id==cpg])
    m=pd.concat(fr); m=m[m.variant_id.isin(gwi.index)].copy()
    if len(m)<10: continue
    m['maf']=np.where(m.af>0.5,1-m.af,m.af)
    # load genotypes for this chrom (cache)
    if chrom not in plink_cache:
        bf=f"/scratch/cy94/rs4477/coloc/ldref/{chrom}"
        plink_cache[chrom]=read_plink1_bin(bf+".bed",bf+".bim",bf+".fam",verbose=False)
    Gb=plink_cache[chrom]
    avail=set(Gb.snp.values)
    snps=[s for s in m.variant_id if s in avail]
    if len(snps)<10: continue
    sub=Gb.sel(snp=snps).values.astype(float)  # samples x snps
    # impute NaN with column mean
    cm=np.nanmean(sub,axis=0); inds=np.where(np.isnan(sub)); sub[inds]=np.take(cm,inds[1])
    ld=np.corrcoef(sub.T)
    m2=m[m.variant_id.isin(snps)].set_index('variant_id').loc[snps].reset_index()
    g=gwi.loc[snps]
    np.savetxt(OUT+f"/{cpg}.ld",ld)
    df=pd.DataFrame({'snp':snps,'mqtl_beta':m2.slope.values,'mqtl_se':m2.slope_se.values,
        'mqtl_maf':m2.maf.values,'mqtl_n':452,'gwas_beta':g.b.values,'gwas_se':g.se.values,
        'gwas_maf':np.where(g.freq.values>0.5,1-g.freq.values,g.freq.values),'gwas_n':g.N.values})
    df.to_csv(OUT+f"/{cpg}.coloc.tsv",sep="\t",index=False); n_ok+=1
    print(cpg,"ok",len(df),flush=True)
print("prepped:",n_ok)
