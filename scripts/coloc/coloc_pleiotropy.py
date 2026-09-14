import pandas as pd, numpy as np, pyarrow.parquet as pq
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from pandas_plink import read_plink1_bin
plt.rcParams.update({"font.size":10,"savefig.dpi":200,"savefig.bbox":"tight","figure.facecolor":"white"})
G="/g/data/cy94/rs4477/downstream/genome_452"
cpg="chr5_56564909"; chrom="chr5"; cpos=56564909; win=200000
# mQTL
pth=f"{G}/{chrom}/nominal/{chrom}.sig.cis_qtl_pairs.{chrom}.parquet"
pf=pq.ParquetFile(pth); fr=[]
for i in range(pf.num_row_groups):
    t=pf.read_row_group(i,columns=['phenotype_id','variant_id','pval_nominal']).to_pandas()
    fr.append(t[t.phenotype_id==cpg])
mq=pd.concat(fr); p=mq.variant_id.str.split(':',expand=True); mq['pos']=p[1].astype(int)
mq=mq[(mq.pos>cpos-win)&(mq.pos<cpos+win)].copy()
# both GWAS
gwas={}
for dis,f in [('CAD','cad'),('T2D','t2d')]:
    g=pd.read_csv(f"/scratch/cy94/rs4477/coloc/gwas/{f}_hg38.ma.full",sep="\t")
    g=g[g.SNP.isin(mq.variant_id)][['SNP','p']].rename(columns={'SNP':'variant_id','p':'gp'})
    gwas[dis]=mq.merge(g,on='variant_id',how='inner')
# LD to lead mQTL
bf=f"/scratch/cy94/rs4477/coloc/ldref/{chrom}"
Gb=read_plink1_bin(bf+".bed",bf+".bim",bf+".fam",verbose=False)
avail=set(Gb.snp.values); mq=mq[mq.variant_id.isin(avail)]
lead=mq.loc[mq.pval_nominal.idxmin(),'variant_id']; snps=list(mq.variant_id)
sub=Gb.sel(snp=snps).values.astype(float); cmn=np.nanmean(sub,axis=0); ii=np.where(np.isnan(sub)); sub[ii]=np.take(cmn,ii[1])
li=snps.index(lead); r2=np.array([np.corrcoef(sub[:,li],sub[:,j])[0,1]**2 for j in range(len(snps))])
mq['r2']=dict(zip(snps,r2)); mq['r2']=[dict(zip(snps,r2))[s] for s in mq.variant_id]
for dis in gwas: gwas[dis]['r2']=gwas[dis].variant_id.map(dict(zip(snps,r2)))
fig,ax=plt.subplots(3,1,figsize=(9,9),sharex=True)
sc=ax[0].scatter(mq.pos/1e6,-np.log10(mq.pval_nominal),c=mq.r2,cmap='RdYlBu_r',vmin=0,vmax=1,s=22,edgecolor='grey',lw=0.2)
ax[0].set_ylabel('-log10 P\n(mQTL)'); ax[0].set_title('chr5:56564909 — pleiotropic methylation locus (CAD + T2D)',fontweight='bold')
for i,dis in enumerate(['CAD','T2D']):
    gg=gwas[dis]
    ax[i+1].scatter(gg.pos/1e6,-np.log10(gg.gp),c=gg.r2,cmap='RdYlBu_r',vmin=0,vmax=1,s=22,edgecolor='grey',lw=0.2)
    ax[i+1].set_ylabel(f'-log10 P\n({dis})')
ax[2].set_xlabel('Chr5 position (Mb)')
plt.colorbar(sc,ax=ax,label='LD (r2 to lead)',fraction=0.03,pad=0.02)
plt.savefig('/scratch/cy94/rs4477/figures_thesis/coloc_pleiotropy_chr5.png'); print('saved',len(mq),'mQTL variants')
