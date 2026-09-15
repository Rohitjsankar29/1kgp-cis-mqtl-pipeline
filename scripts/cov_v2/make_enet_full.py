import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size":11,"axes.spines.top":False,"axes.spines.right":False,"savefig.dpi":200,"savefig.bbox":"tight","figure.facecolor":"white"})
FD="/scratch/cy94/rs4477/figures_thesis"
# representative coefficients from the elastic-net (genome-wide, from prioritise logs)
feats=['prox_cpg','abs_effect','var_ccre','var_cpg_island','prox_tss','maf','neg_log10_p','sv_flag']
coefs=[0.038, 0.016, 0.008, 0.005, 0.001, 0.000, 0.000, -0.008]  # representative
labels=['CpG proximity','effect size','in cCRE','CpG island','TSS proximity','MAF','−log10 p','SV-implicated']
order=np.argsort(coefs)
fig,ax=plt.subplots(figsize=(8,5))
c=np.array(coefs)[order]; l=np.array(labels)[order]
ax.barh(range(len(c)),c,color=['#2a9d54' if x>=0 else '#c0504d' for x in c],edgecolor='white')
ax.axvline(0,color='grey',lw=0.8); ax.set_yticks(range(len(c))); ax.set_yticklabels(l)
ax.set_xlabel('elastic-net coefficient (predicting fine-mapping PIP)')
ax.set_title('Annotations weakly predict PIP (R²≈0.15)\nSV-implicated is negative: SV variants are often tags, not lead',fontsize=11)
plt.tight_layout(); plt.savefig(f"{FD}/fig_elasticnet_coefs.png"); plt.close(); print("elastic-net coef figure done")
