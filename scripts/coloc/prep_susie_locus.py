import sys, pandas as pd, numpy as np, os, pyarrow.parquet as pq, subprocess
G="/g/data/cy94/rs4477/downstream/genome_452"
OUT="/scratch/cy94/rs4477/coloc/susie_loci"; os.makedirs(OUT,exist_ok=True)
PLINK="/scratch/cy94/rs4477/micromamba/envs/plink/bin/plink2"
# read SMR+HEIDI hits
hits=pd.read_csv("/scratch/cy94/rs4477/coloc/eGFR_smr_heidi_hits.tsv",sep="\t")
# eGFR GWAS (hg38 lifted)
gw=pd.read_csv("/scratch/cy94/rs4477/coloc/gwas/eGFR_hg38.ma.full",sep="\t")
gw['SNP']=gw.SNP.astype(str)
gwi=gw.set_index('SNP')
for _,hit in hits.iterrows():
    cpg=hit.probeID; chrom=cpg.split('_')[0]; chrn=chrom.replace('chr','')
    pth=f"{G}/{chrom}/nominal/{chrom}.sig.cis_qtl_pairs.{chrom}.parquet"
    if chrom in ['chr1','chr2']: pth=f"{G}/{chrom}/nominal/{chrom}.sig.merged.parquet"
    pf=pq.ParquetFile(pth); frames=[]
    for i in range(pf.num_row_groups):
        t=pf.read_row_group(i,columns=['phenotype_id','variant_id','af','slope','slope_se']).to_pandas()
        frames.append(t[t.phenotype_id==cpg])
    m=pd.concat(frames)
    if len(m)<10: continue
    m['maf']=np.where(m.af>0.5,1-m.af,m.af)
    # match to GWAS by variant_id (both hg38 chr:pos:ref:alt); GWAS has snp_fwd/rev
    m['g_fwd']=m.variant_id.isin(gwi.index)
    mm=m[m.g_fwd].copy()
    if len(mm)<10: continue
    g=gwi.loc[mm.variant_id]
    # LD matrix via plink (r) for these variants
    snplist=OUT+f"/{cpg}.snps"; mm.variant_id.to_csv(snplist,index=False,header=False)
    bfile=f"/scratch/cy94/rs4477/coloc/ldref/{chrom}"
    subprocess.run([PLINK,"--bfile",bfile,"--extract",snplist,"--r-unphased","square","--out",OUT+f"/{cpg}"],capture_output=True)
    import glob
    lf=glob.glob(OUT+f"/{cpg}*.vcor1")+glob.glob(OUT+f"/{cpg}.ld")+glob.glob(OUT+f"/{cpg}*.rel")
    ldf=lf[0] if lf else OUT+f"/{cpg}.ld"
    if not os.path.exists(ldf): continue
    # coloc input (aligned)
    df=pd.DataFrame({'snp':mm.variant_id.values,'mqtl_beta':mm.slope.values,'mqtl_se':mm.slope_se.values,
        'mqtl_maf':mm.maf.values,'mqtl_n':452,
        'gwas_beta':g.b.values,'gwas_se':g.se.values,
        'gwas_maf':np.where(g.freq.values>0.5,1-g.freq.values,g.freq.values),'gwas_n':g.N.values})
    df.to_csv(OUT+f"/{cpg}.coloc.tsv",sep="\t",index=False)
    print(cpg, "prepped", len(df), flush=True)
print("DONE")
