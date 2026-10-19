import gzip, glob, math
BESD_GLOB="/scratch/cy94/rs4477/coloc/besd/mqtl_chr*.epi"
PRIOR_GZ ="/g/data/cy94/rs4477/downstream/genome_452/genome.prioritised_final_asm_rescored.txt.gz"
HIT_FILES={"eGFR":"/scratch/cy94/rs4477/coloc/eGFR_smr_heidi_hits.tsv",
           "CAD":"/scratch/cy94/rs4477/coloc/cvd_smr_heidi_hits.tsv",
           "T2D":"/scratch/cy94/rs4477/coloc/t2d_smr_heidi_hits.tsv"}
CPGCOL,SCCOL=2,37
def lchoose(n,k):
    if k<0 or k>n: return float("-inf")
    return math.lgamma(n+1)-math.lgamma(k+1)-math.lgamma(n-k+1)
def hsf(k,N,K,n):
    lo,hi=max(0,n-(N-K)),min(n,K); k=max(k,lo); den=lchoose(N,n); t=0.0
    for i in range(k,hi+1): t+=math.exp(lchoose(K,i)+lchoose(N-K,n-i)-den)
    return min(1.0,t)
def nsf2(z): return math.erfc(abs(z)/math.sqrt(2))
besd=set()
for f in glob.glob(BESD_GLOB):
    for l in open(f):
        p=l.split()
        if len(p)>1: besd.add(p[1])
print("BESD CpGs:",len(besd))
mx={}
with gzip.open(PRIOR_GZ,"rt") as fh:
    fh.readline()
    for l in fh:
        c=l.rstrip("\n").split("\t")
        if c[CPGCOL-1] in besd:
            s=int(float(c[SCCOL-1]))
            if s>mx.get(c[CPGCOL-1],-1): mx[c[CPGCOL-1]]=s
hits=set()
for path in HIT_FILES.values():
    fh=open(path); fh.readline()
    for l in fh:
        cpg=l.split("\t")[0]
        if cpg: hits.add(cpg)
print("disease-hit CpGs:",len(hits))
scores=sorted({s for s in mx.values()})
bg ={s:sum(1 for v in mx.values() if v==s) for s in scores}
dis={s:sum(1 for c in hits if mx.get(c)==s) for s in scores}
N=sum(bg.values()); n=sum(dis.values()); overall=n/N
with open("enrichment_by_convergence_stats.tsv","w") as o:
    o.write("score\tbg\tdisease\trate_pct\tfold\n")
    for s in scores:
        r=dis[s]/bg[s] if bg[s] else 0
        o.write(f"{s}\t{bg[s]}\t{dis[s]}\t{r*100:.4f}\t{r/overall:.3f}\n")
    for thr in sorted(scores,reverse=True):
        K=sum(bg[s] for s in scores if s>=thr); k=sum(dis[s] for s in scores if s>=thr)
        o.write(f"# >= {thr}: {k}/{n} vs {K}/{N} fold={(k/n)/(K/N):.2f} P={hsf(k,N,K,n):.3e}\n")
    T=sum(dis[s]*s for s in scores)-n*sum(bg[s]*s for s in scores)/N
    var=overall*(1-overall)*(sum(bg[s]*s*s for s in scores)-(sum(bg[s]*s for s in scores))**2/N)
    z=T/math.sqrt(var)
    o.write(f"# CA trend z={z:.2f} P={nsf2(z):.3e}\n")
print("CA trend z=%.2f P=%.2e"%(z,nsf2(z)))
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
rate=[dis[s]/bg[s]*100 for s in scores]
ramp=["#bcdcec","#7bb8da","#357fb0","#0f4c75","#07314d"]
cols=ramp[-len(scores):]
fig,ax=plt.subplots(figsize=(5.4,4.1),dpi=300)
ax.bar(scores,rate,width=0.72,color=cols,edgecolor="white",linewidth=1.2,zorder=3)
ax.axhline(overall*100,ls=(0,(4,3)),lw=1.2,color="#9aa0a6",zorder=2)
ax.text(max(scores)+0.46,overall*100,"all tested\n%.2f%%"%(overall*100),va="center",ha="left",fontsize=7.5,color="#6b7177")
for s,r in zip(scores,rate):
    ax.text(s,r+0.012,"%.2f%%"%r,ha="center",va="bottom",fontsize=9,fontweight="bold",color="#223")
    ax.text(s,-0.028,"%d/%s"%(dis[s],format(bg[s],",")),ha="center",va="top",fontsize=7,color="#6b7177")
ax.text(0.03,0.95,"Cochran\u2013Armitage trend\n$z$ = %.1f,  $P$ = %.0e"%(z,nsf2(z)),transform=ax.transAxes,va="top",ha="left",fontsize=8.5,bbox=dict(boxstyle="round,pad=0.4",fc="#f4f6f8",ec="#d7dbdf",lw=0.8))
ax.set_xlabel("Convergence score (evidence axes met)",fontsize=10)
ax.set_ylabel("cis-mQTLs mediating disease (%)",fontsize=10)
ax.set_title("Disease-mediating loci are enriched at higher convergence",fontsize=10.5,fontweight="bold",pad=10)
ax.set_xticks(scores); ax.set_ylim(0,max(rate)*1.18); ax.set_xlim(min(scores)-0.6,max(scores)+0.9)
ax.yaxis.set_major_formatter(FuncFormatter(lambda v,_:"%.2f"%v))
for sp in ["top","right"]: ax.spines[sp].set_visible(False)
ax.spines["left"].set_color("#c7ccd1"); ax.spines["bottom"].set_color("#c7ccd1")
ax.tick_params(colors="#55595e",labelsize=8.5); ax.grid(axis="y",lw=0.5,color="#eceef0",zorder=0)
plt.tight_layout(); plt.savefig("fig_disease_enrichment_by_convergence.png",bbox_inches="tight",facecolor="white")
print("wrote fig_disease_enrichment_by_convergence.png")
