import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size":11,"axes.spines.top":False,"axes.spines.right":False,
                     "savefig.dpi":200,"savefig.bbox":"tight","figure.facecolor":"white"})
G="/g/data/cy94/rs4477/downstream/genome_452"; FD="/scratch/cy94/rs4477/figures_thesis"

# === FIG 3: Co-methylation domains ===
v=pd.read_csv(G+"/genome.prioritised_final_asm_classified.txt.gz",sep="\t",
              usecols=lambda c:c in['n_cpgs','sv_implicated','asm_frac_fisher'])
fig,ax=plt.subplots(1,3,figsize=(15,4.2))
ax[0].hist(np.log10(v.n_cpgs),bins=50,color="#2c6fbb",edgecolor="white")
ax[0].set_xlabel("log10(domain size, CpGs)"); ax[0].set_ylabel("variants"); ax[0].set_title("Co-methylation domain sizes")
sv=v[v.sv_implicated==True].n_cpgs; nosv=v[v.sv_implicated==False].n_cpgs
bp=ax[1].boxplot([np.log10(nosv),np.log10(sv)],tick_labels=["no SV","SV"],showfliers=False,patch_artist=True)
bp['boxes'][0].set_facecolor("#9aa0a6"); bp['boxes'][1].set_facecolor("#e07b39")
ax[1].set_ylabel("log10(domain size)"); ax[1].set_title("SV variants drive larger domains\n(9.7 vs 3.9 CpGs, p<1e-300)")
vt=v[v.asm_frac_fisher.notna()].copy()
vt['bin']=pd.cut(vt.n_cpgs,[0,5,20,50,100,10000],labels=['1-5','6-20','21-50','51-100','100+'])
med=vt.groupby('bin',observed=True).asm_frac_fisher.mean()
ax[2].bar(range(len(med)),med.values,color="#2a9d54")
ax[2].set_xticks(range(len(med))); ax[2].set_xticklabels(med.index)
ax[2].set_xlabel("domain size (CpGs)"); ax[2].set_ylabel("mean ASM fraction"); ax[2].set_title("Larger domains show more ASM\n(rho=0.19)")
plt.tight_layout(); plt.savefig(f"{FD}/fig3_comethyl_domains.png"); plt.close(); print("fig3 domains done")

# === FIG 4: SMR cardiometabolic (per-disease, top hits) ===
diseases={'eGFR':'eGFR (kidney)','cvd':'CAD (coronary)','t2d':'T2D (diabetes)'}
fig,ax=plt.subplots(1,3,figsize=(16,4.5))
for i,(d,lab) in enumerate(diseases.items()):
    h=pd.read_csv(f"/scratch/cy94/rs4477/coloc/{d}_smr_heidi_hits.tsv",sep="\t")
    h=h.sort_values('p_SMR').head(12)
    genes=[g if not str(g).startswith('ENSG') else p.split('_')[0]+':'+p.split('_')[1][:6] for g,p in zip(h.nearest_gene,h.probeID)]
    y=range(len(h))
    ax[i].barh(y,-np.log10(h.p_SMR),color='#c0504d',edgecolor='white')
    ax[i].set_yticks(y); ax[i].set_yticklabels(genes,fontsize=8); ax[i].invert_yaxis()
    ax[i].set_xlabel('-log10(SMR p)'); ax[i].set_title(f'{lab}\n{len(pd.read_csv(f"/scratch/cy94/rs4477/coloc/{d}_smr_heidi_hits.tsv",sep=chr(9)))} mediated loci')
plt.tight_layout(); plt.savefig(f"{FD}/fig4_smr_cardiometabolic.png"); plt.close(); print("fig4 SMR done")
