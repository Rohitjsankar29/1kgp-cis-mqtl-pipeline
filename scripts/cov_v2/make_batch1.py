import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size":11,"axes.spines.top":False,"axes.spines.right":False,"savefig.dpi":200,"savefig.bbox":"tight","figure.facecolor":"white"})
G="/g/data/cy94/rs4477/downstream/genome_452"; FD="/scratch/cy94/rs4477/figures_thesis"
v=pd.read_csv(G+"/genome.prioritised_final_asm_classified.txt.gz",sep="\t",
              usecols=['pip','n_cpgs','ccre','sv_implicated','score_convergence','asm_frac_fisher','asm_abs_delta','imprinted','nearest_gene','variant_id'])

# FIG A: PIP distribution
fig,ax=plt.subplots(figsize=(7,4.5))
ax.hist(v.pip,bins=50,color='#2c6fbb',edgecolor='white')
ax.axvline(0.9,color='red',ls='--',label='PIP>0.9 (credible-causal)')
ax.set_xlabel('fine-mapping PIP'); ax.set_ylabel('variants'); ax.set_yscale('log')
ax.set_title('Fine-mapping posterior inclusion probabilities'); ax.legend()
plt.savefig(f"{FD}/fig_pip_dist.png"); plt.close(); print("PIP dist done")

# FIG B: Functional enrichment (cCRE classes)
cc=v.ccre.astype(str).value_counts()
cc=cc[cc.index!='none'].head(8)
fig,ax=plt.subplots(figsize=(7,4.5))
ax.barh(range(len(cc)),cc.values,color='#2a9d54',edgecolor='white')
ax.set_yticks(range(len(cc))); ax.set_yticklabels(cc.index); ax.invert_yaxis()
ax.set_xlabel('prioritised variants'); ax.set_title('Regulatory context: cCRE class enrichment')
plt.savefig(f"{FD}/fig_functional_enrich.png"); plt.close(); print("functional done")

# FIG C: Convergence score distribution
sc=v.score_convergence.value_counts().sort_index()
fig,ax=plt.subplots(figsize=(7,4.5))
ax.bar(sc.index,sc.values,color=['#cccccc','#9ec5e8','#6ba3d6','#2c6fbb','#e07b39','#c0504d'][:len(sc)],edgecolor='white')
ax.set_xlabel('convergence score (independent evidence axes)'); ax.set_ylabel('variants'); ax.set_yscale('log')
ax.set_title('Multi-evidence convergence: variants by number of supporting axes')
for i,val in zip(sc.index,sc.values): ax.text(i,val*1.1,f'{val:,}',ha='center',fontsize=8)
plt.savefig(f"{FD}/fig_convergence.png"); plt.close(); print("convergence done")

# FIG D: n_cpgs (co-methylation domain) + effect
fig,ax=plt.subplots(1,2,figsize=(11,4.2))
ax[0].hist(np.log10(v.n_cpgs),bins=50,color='#2c6fbb',edgecolor='white')
ax[0].set_xlabel('log10(CpGs per variant)'); ax[0].set_ylabel('variants'); ax[0].set_title('CpGs controlled per variant')
# PIP vs ASM colored (the orthogonality, cleaner)
vt=v[v.asm_frac_fisher.notna()].sample(min(8000,v.asm_frac_fisher.notna().sum()),random_state=1)
ax[1].scatter(vt.pip,vt.asm_frac_fisher,s=4,alpha=0.2,color='#8B2252')
ax[1].set_xlabel('PIP'); ax[1].set_ylabel('ASM fraction'); ax[1].set_title('ASM vs PIP (r=0.011)')
plt.tight_layout(); plt.savefig(f"{FD}/fig_ncpg_asm.png"); plt.close(); print("ncpg/asm done")

# FIG E: Imprinting recovery (imprinted vs sequence-driven ASM)
imp=v[(v.imprinted==True)&(v.asm_frac_fisher.notna())]
seq=v[(v.imprinted==False)&(v.asm_frac_fisher.notna())]
fig,ax=plt.subplots(figsize=(7,4.5))
ax.hist(seq.asm_frac_fisher,bins=40,alpha=0.6,label=f'sequence-driven (n={len(seq)})',color='#2c6fbb',density=True)
ax.hist(imp.asm_frac_fisher,bins=40,alpha=0.6,label=f'imprinted control (n={len(imp)})',color='#c0504d',density=True)
ax.set_xlabel('ASM fraction (Fisher)'); ax.set_ylabel('density')
ax.set_title('ASM: imprinted loci (positive control) vs sequence-driven'); ax.legend()
plt.savefig(f"{FD}/fig_imprinting.png"); plt.close(); print("imprinting done")
print("BATCH 1 COMPLETE")
