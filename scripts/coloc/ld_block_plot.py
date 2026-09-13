import sys, pandas as pd, numpy as np, os, pyarrow.parquet as pq
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from pandas_plink import read_plink1_bin
plt.rcParams.update({"font.size":10,"savefig.dpi":200,"savefig.bbox":"tight","figure.facecolor":"white"})
G="/g/data/cy94/rs4477/downstream/genome_452"
focal=sys.argv[1]; cpg=sys.argv[2]; gene=sys.argv[3]; out=sys.argv[4]; win=int(sys.argv[5]) if len(sys.argv)>5 else 50000
chrom=focal.split(':')[0]; fpos=int(focal.split(':')[1])
pth=f"{G}/{chrom}/nominal/{chrom}.sig.cis_qtl_pairs.{chrom}.parquet"
if chrom in ['chr1','chr2']: pth=f"{G}/{chrom}/nominal/{chrom}.sig.merged.parquet"
pf=pq.ParquetFile(pth); fr=[]
for i in range(pf.num_row_groups):
    t=pf.read_row_group(i,columns=['phenotype_id','variant_id','pval_nominal','slope']).to_pandas()
    fr.append(t[t.phenotype_id==cpg])
m=pd.concat(fr)
p=m.variant_id.str.split(':',expand=True); m['pos']=p[1].astype(int)
m=m[(m.pos>fpos-win)&(m.pos<fpos+win)].sort_values('pos').reset_index(drop=True)
if len(m)<3: print("too few variants:",len(m)); sys.exit()
bf=f"/scratch/cy94/rs4477/coloc/ldref/{chrom}"
Gb=read_plink1_bin(bf+".bed",bf+".bim",bf+".fam",verbose=False)
avail=set(Gb.snp.values); m=m[m.variant_id.isin(avail)].reset_index(drop=True)
if len(m)<3: print("too few after LD filter"); sys.exit()
sub=Gb.sel(snp=list(m.variant_id)).values.astype(float)
cmn=np.nanmean(sub,axis=0); inds=np.where(np.isnan(sub)); sub[inds]=np.take(cmn,inds[1])
ld=np.corrcoef(sub.T)**2
cl=pd.read_csv(G+"/genome.prioritised_final_asm_classified.txt.gz",sep="\t",usecols=['variant_id','pip','asm_frac_fisher'])
m=m.merge(cl,on='variant_id',how='left')
n=len(m); x=np.arange(n)
fig=plt.figure(figsize=(11,9)); gs=GridSpec(4,1,height_ratios=[1.2,1,1,4],hspace=0.15)
fi=list(m.variant_id).index(focal) if focal in list(m.variant_id) else int(np.argmax(-np.log10(m.pval_nominal)))
ax0=fig.add_subplot(gs[0]); ax0.scatter(x,-np.log10(m.pval_nominal),c='#8B2252',s=15)
ax0.scatter(fi,-np.log10(m.pval_nominal.iloc[fi]),c='red',s=70,marker='D',zorder=5)
ax0.set_ylabel('-log10 P\n(mQTL)'); ax0.set_xticks([]); ax0.set_title(f'{gene} locus  (focal {focal}, PIP={m.pip.iloc[fi]:.2f}, ASM={m.asm_frac_fisher.iloc[fi]:.2f})',fontweight='bold',fontsize=11)
ax1=fig.add_subplot(gs[1]); ax1.bar(x,m.pip.fillna(0),color='#2c6fbb'); ax1.set_ylabel('PIP'); ax1.set_xticks([]); ax1.set_ylim(0,1)
ax2=fig.add_subplot(gs[2]); ax2.bar(x,m.asm_frac_fisher.fillna(0),color='#2a9d54'); ax2.set_ylabel('ASM\nfraction'); ax2.set_xticks([]); ax2.set_ylim(0,1)
ax3=fig.add_subplot(gs[3]); im=ax3.imshow(ld,cmap='Reds',aspect='auto',vmin=0,vmax=1)
ax3.axvline(fi,color='blue',lw=0.8); ax3.axhline(fi,color='blue',lw=0.8)
ax3.set_xlabel('variants (ordered by position)'); ax3.set_ylabel('variants')
plt.colorbar(im,ax=ax3,label='LD (r2)',fraction=0.046)
plt.savefig(out); print(f"saved {out}: {n} variants, focal PIP={m.pip.iloc[fi]:.2f} ASM={m.asm_frac_fisher.iloc[fi]:.2f}")
