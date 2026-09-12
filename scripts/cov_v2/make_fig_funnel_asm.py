import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size":11,"axes.spines.top":False,"axes.spines.right":False,
                     "savefig.dpi":200,"savefig.bbox":"tight","figure.facecolor":"white"})
G="/g/data/cy94/rs4477/downstream/genome_452"
FD="/scratch/cy94/rs4477/figures_thesis"

# === FIG 1: Discovery funnel ===
stages=['CpGs\ntested','Significant\nmQTLs\n(FDR<0.05)','Credible-set\nvariants','Credible\ncausal\n(PIP>0.9)']
vals=[25912587, 1638801, 4254269, 192172]
fig,ax=plt.subplots(figsize=(7,4.5))
colors=['#2c6fbb','#3a8fb7','#5aa9c9','#7cc4d6']
bars=ax.bar(range(len(stages)),vals,color=colors,edgecolor='white',width=0.65)
ax.set_yscale('log'); ax.set_xticks(range(len(stages))); ax.set_xticklabels(stages,fontsize=9)
ax.set_ylabel('count (log scale)'); ax.set_title('Genome-wide cis-mQTL discovery funnel')
for b,v in zip(bars,vals): ax.text(b.get_x()+b.get_width()/2,v*1.15,f'{v:,}',ha='center',fontsize=8)
plt.savefig(f"{FD}/fig1_discovery_funnel.png"); plt.close()
print("fig1 funnel done")

# === FIG 2: ASM vs PIP (the orthogonality finding) ===
# stratified data
try:
    st=pd.read_csv(G+"/genome.asm_by_stratum.tsv",sep="\t")
    order=["PIP>0.9","PIP_0.5-0.9","PIP_0.1-0.5","PIP<0.1"]
    fig,ax=plt.subplots(1,2,figsize=(11,4.2))
    # panel A: mean ASM fraction by stratum
    means=[st[st.stratum==s].frac_mwu_sig.mean() for s in order]
    ax[0].bar(range(4),means,color='#e07b39',edgecolor='white')
    ax[0].set_xticks(range(4)); ax[0].set_xticklabels(['>0.9','0.5-0.9','0.1-0.5','<0.1'])
    ax[0].set_xlabel('fine-mapping PIP stratum'); ax[0].set_ylabel('mean ASM fraction (MWU)')
    ax[0].set_title('ASM is flat across PIP strata\n(matched on mQTL strength)')
    # panel B: scatter ASM vs PIP (from classified table)
    v=pd.read_csv(G+"/genome.prioritised_final_asm_classified.txt.gz",sep="\t",usecols=lambda c:c in['pip','asm_frac_fisher'])
    v=v.dropna(subset=['asm_frac_fisher']).sample(min(5000,len(v)),random_state=1)
    ax[1].scatter(v.pip,v.asm_frac_fisher,s=3,alpha=0.15,color='#2c6fbb')
    ax[1].set_xlabel('fine-mapping PIP'); ax[1].set_ylabel('ASM fraction (Fisher)')
    ax[1].set_title('ASM uncorrelated with PIP\n(r=0.011)')
    plt.tight_layout(); plt.savefig(f"{FD}/fig2_asm_pip.png"); plt.close()
    print("fig2 ASM-PIP done")
except Exception as e:
    print("fig2 error:",e)
