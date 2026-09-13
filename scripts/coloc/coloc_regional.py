import sys, pandas as pd, numpy as np, os, pyarrow.parquet as pq
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from pandas_plink import read_plink1_bin
plt.rcParams.update({"font.size":10,"savefig.dpi":200,"savefig.bbox":"tight","figure.facecolor":"white"})
G="/g/data/cy94/rs4477/downstream/genome_452"
cpg=sys.argv[1]; gene=sys.argv[2]; disease=sys.argv[3]; out=sys.argv[4]; win=int(sys.argv[5]) if len(sys.argv)>5 else 200000
chrom=cpg.split('_')[0]; cpos=int(cpg.split('_')[1]); chrn=chrom.replace('chr','')
# mQTL association for this CpG
pth=f"{G}/{chrom}/nominal/{chrom}.sig.cis_qtl_pairs.{chrom}.parquet"
if chrom in ['chr1','chr2']: pth=f"{G}/{chrom}/nominal/{chrom}.sig.merged.parquet"
pf=pq.ParquetFile(pth); fr=[]
for i in range(pf.num_row_groups):
    t=pf.read_row_group(i,columns=['phenotype_id','variant_id','pval_nominal']).to_pandas()
    fr.append(t[t.phenotype_id==cpg])
mq=pd.concat(fr)
p=mq.variant_id.str.split(':',expand=True); mq['pos']=p[1].astype(int)
mq=mq[(mq.pos>cpos-win)&(mq.pos<cpos+win)].copy()
# GWAS for this region
gwf={'eGFR':'eGFR','cvd':'cad','t2d':'t2d'}
gw=pd.read_csv(f"/scratch/cy94/rs4477/coloc/gwas/{gwf[disease]}_hg38.ma.full",sep="\t")
gw=gw[gw.SNP.isin(mq.variant_id)][['SNP','p']].rename(columns={'SNP':'variant_id','p':'gwas_p'})
mq=mq.merge(gw,on='variant_id',how='inner')
if len(mq)<5: print("too few matched"); sys.exit()
# LD to lead mQTL variant
bf=f"/scratch/cy94/rs4477/coloc/ldref/{chrom}"
Gb=read_plink1_bin(bf+".bed",bf+".bim",bf+".fam",verbose=False)
avail=set(Gb.snp.values); mq=mq[mq.variant_id.isin(avail)].copy()
lead=mq.loc[mq.pval_nominal.idxmin(),'variant_id']
snps=list(mq.variant_id)
sub=Gb.sel(snp=snps).values.astype(float)
cmn=np.nanmean(sub,axis=0); inds=np.where(np.isnan(sub)); sub[inds]=np.take(cmn,inds[1])
li=snps.index(lead)
r2=np.array([np.corrcoef(sub[:,li],sub[:,j])[0,1]**2 for j in range(len(snps))])
mq['r2']=r2
# two-panel: mQTL (top) + GWAS (bottom), colored by LD
fig,ax=plt.subplots(2,1,figsize=(10,7),sharex=True)
sc=ax[0].scatter(mq.pos/1e6,-np.log10(mq.pval_nominal),c=mq.r2,cmap='RdYlBu_r',vmin=0,vmax=1,s=25,edgecolor='grey',linewidth=0.2)
ax[0].scatter(mq.loc[mq.variant_id==lead,'pos']/1e6,-np.log10(mq.pval_nominal.min()),marker='D',c='purple',s=80,zorder=5,label='lead')
ax[0].set_ylabel('-log10 P (mQTL)'); ax[0].set_title(f'{gene}: methylation QTL vs {disease.upper()} colocalization',fontweight='bold')
ax[0].legend()
ax[1].scatter(mq.pos/1e6,-np.log10(mq.gwas_p),c=mq.r2,cmap='RdYlBu_r',vmin=0,vmax=1,s=25,edgecolor='grey',linewidth=0.2)
ax[1].set_ylabel(f'-log10 P ({disease.upper()})'); ax[1].set_xlabel(f'Chr{chrn} position (Mb)')
plt.colorbar(sc,ax=ax,label='LD (r2 to lead)',fraction=0.03,pad=0.02)
plt.savefig(out); print(f"saved {out}: {len(mq)} variants, lead={lead}")
