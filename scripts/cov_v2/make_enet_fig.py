import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from sklearn.linear_model import ElasticNetCV
plt.rcParams.update({"font.size":11,"axes.spines.top":False,"axes.spines.right":False,"savefig.dpi":200,"savefig.bbox":"tight","figure.facecolor":"white"})
G="/g/data/cy94/rs4477/downstream/genome_452"; FD="/scratch/cy94/rs4477/figures_thesis"
v=pd.read_csv(G+"/genome.prioritised_final_asm_classified.txt.gz",sep="\t",
              usecols=['pip','n_cpgs','ccre','sv_implicated','asm_frac_fisher','asm_abs_delta'])
v['in_ccre']=(v.ccre.astype(str)!='none').astype(float)
v['sv']=v.sv_implicated.fillna(False).astype(float)
v['log_ncpg']=np.log1p(v.n_cpgs)
feats1=['in_ccre','sv','log_ncpg']
feats_asm=['in_ccre','sv','log_ncpg','asm_frac_fisher','asm_abs_delta']
X1=v[feats1].fillna(0).values; y=v['pip'].values
en1=ElasticNetCV(cv=5,max_iter=3000,n_jobs=-1).fit(X1,y); r2_1=en1.score(X1,y)
vt=v[v.asm_frac_fisher.notna()].copy()
X2=vt[feats_asm].fillna(0).values; y2=vt['pip'].values
en2=ElasticNetCV(cv=5,max_iter=3000,n_jobs=-1).fit(X2,y2); r2_2=en2.score(X2,y2)
print(f"Model1 R2={r2_1:.3f}", dict(zip(feats1,np.round(en1.coef_,4))))
print(f"Model2 R2={r2_2:.3f}", dict(zip(feats_asm,np.round(en2.coef_,4))))
fig,ax=plt.subplots(1,2,figsize=(12,4.5))
n1=['in cCRE','SV-impl.','log(n CpGs)']; c1=en1.coef_
ax[0].barh(range(len(c1)),c1,color=['#2a9d54' if x>=0 else '#c0504d' for x in c1],edgecolor='white')
ax[0].axvline(0,color='grey',lw=0.8); ax[0].set_yticks(range(len(c1))); ax[0].set_yticklabels(n1)
ax[0].set_xlabel('elastic-net coefficient'); ax[0].set_title(f'Annotations predict PIP poorly (R2={r2_1:.2f})')
n2=['in cCRE','SV-impl.','log(n CpGs)','ASM fraction','ASM |delta|']; c2=en2.coef_
ax[1].barh(range(len(c2)),c2,color=['#2a9d54' if x>=0 else '#c0504d' for x in c2],edgecolor='white')
ax[1].axvline(0,color='grey',lw=0.8); ax[1].set_yticks(range(len(c2))); ax[1].set_yticklabels(n2)
ax[1].set_xlabel('elastic-net coefficient'); ax[1].set_title(f'Adding ASM does not help (R2={r2_2:.2f})')
plt.tight_layout(); plt.savefig(f"{FD}/fig_elasticnet.png"); plt.close(); print("done")
