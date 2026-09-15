import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size":10,"axes.spines.top":False,"axes.spines.right":False,"savefig.dpi":200,"savefig.bbox":"tight","figure.facecolor":"white"})
FD="/scratch/cy94/rs4477/figures_thesis"
# SMR forest: top mediated loci per disease with b_SMR +/- CI
fig,ax=plt.subplots(1,3,figsize=(15,5))
diseases={'eGFR':('eGFR (kidney)','#2c6fbb'),'cvd':('CAD (coronary)','#c0504d'),'t2d':('T2D (diabetes)','#e07b39')}
for i,(d,(lab,col)) in enumerate(diseases.items()):
    h=pd.read_csv(f'/scratch/cy94/rs4477/coloc/{d}_smr_heidi_hits.tsv',sep='\t')
    h=h.sort_values('p_SMR').head(12).iloc[::-1]
    genes=[g if not str(g).startswith('ENSG') else p.split('_')[0]+':'+p.split('_')[1][:6] for g,p in zip(h.nearest_gene,h.probeID)]
    y=range(len(h))
    ax[i].errorbar(h.b_SMR, y, xerr=1.96*h.se_SMR, fmt='o', color=col, capsize=3, markersize=5)
    ax[i].axvline(0,color='grey',ls='--',lw=0.8)
    ax[i].set_yticks(y); ax[i].set_yticklabels(genes,fontsize=8)
    ax[i].set_xlabel('SMR effect (b_SMR) ± 95% CI'); ax[i].set_title(f'{lab}\n{len(pd.read_csv(f"/scratch/cy94/rs4477/coloc/{d}_smr_heidi_hits.tsv",sep=chr(9)))} mediated loci')
plt.tight_layout(); plt.savefig(f"{FD}/fig_smr_forest.png"); plt.close(); print("SMR forest done")
