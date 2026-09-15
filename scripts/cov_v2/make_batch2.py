import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size":11,"axes.spines.top":False,"axes.spines.right":False,"savefig.dpi":200,"savefig.bbox":"tight","figure.facecolor":"white"})
FD="/scratch/cy94/rs4477/figures_thesis"
cov=pd.read_csv("/scratch/cy94/rs4477/1kgp-cis-mqtl/covariates/genome_452/chr20.covariates.tsv",sep="\t")
covm=cov.set_index("ID").T; covm.index.name="sample"
# 1000G panel: sample pop super_pop gender
panel=pd.read_csv("/scratch/cy94/rs4477/reference/integrated_call_samples_v3.20130502.ALL.panel",sep="\t")
pc=panel.columns
# figure out the superpop column
spcol=[c for c in pc if 'super' in c.lower()]
scol=spcol[0] if spcol else pc[2]
idcol=pc[0]
popmap=dict(zip(panel[idcol],panel[scol]))
# sample IDs may be GM-prefixed vs HG - strip GM
def norm(s):
    return s.replace('GM','NA') if s.startswith('GM') else s
covm['pop']=[popmap.get(s, popmap.get(norm(s),'other')) for s in covm.index]
print('population counts:'); print(covm['pop'].value_counts())
fig,ax=plt.subplots(figsize=(7.5,6))
colors={'EUR':'#2c6fbb','AFR':'#e07b39','EAS':'#2a9d54','SAS':'#c0504d','AMR':'#7b5ea7','other':'#999999'}
for pop in covm['pop'].value_counts().index:
    sub=covm[covm['pop']==pop]
    ax.scatter(sub['genoPC1'].astype(float),sub['genoPC2'].astype(float),s=28,label=f'{pop} (n={len(sub)})',
               alpha=0.75,color=colors.get(pop,'#999999'),edgecolor='white',linewidth=0.3)
ax.set_xlabel('Genotype PC1'); ax.set_ylabel('Genotype PC2')
ax.set_title(f'Multi-ancestry cohort structure (n={len(covm)})'); ax.legend(title='superpopulation',fontsize=9)
plt.savefig(f"{FD}/fig_ancestry_pca.png"); plt.close(); print("ancestry PCA done")
