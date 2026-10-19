import pandas as pd, numpy as np, os
from pyliftover import LiftOver
import pyarrow.parquet as pq
lo=LiftOver('hg38','hg19')
G="/g/data/cy94/rs4477/downstream/genome_452"

def coloc_abf(b1,v1,b2,v2, p1=1e-4,p2=1e-4,p12=1e-5):
    # per-SNP log ABF for each trait (Wakefield), then combine
    def labf(b,v,W=0.15**2):
        r=W/(W+v); z2=b*b/v
        return 0.5*(np.log(1-r))+0.5*r*z2
    l1=labf(b1,v1); l2=labf(b2,v2)
    # log sum exp helper
    def lse(x):
        m=np.max(x); return m+np.log(np.sum(np.exp(x-m)))
    n=len(b1)
    lH1=np.log(p1)+lse(l1)
    lH2=np.log(p2)+lse(l2)
    lH3=np.log(p1)+np.log(p2)+lse(l1)+lse(l2)  # approx (not excluding same-snp)
    # H3 proper: sum over i!=j -> lse(l1)+lse(l2) minus same-index; approximate
    lH4=np.log(p12)+lse(l1+l2)
    lH0=0.0
    # normalize
    arr=np.array([lH0,lH1,lH2,lH3,lH4])
    m=np.max(arr); post=np.exp(arr-m); post/=post.sum()
    return post  # [H0,H1,H2,H3,H4]

# load eGFR once
print("loading eGFR...",flush=True)
gw=pd.read_csv("gwas/eGFR_stanzick2021.gz",sep='\t',compression='gzip',usecols=['chr','pos','Allele1','Allele2','Freq1','Effect','StdErr'])
gw['key']=gw.chr.astype(str)+':'+gw.pos.astype(str)
gw=gw.drop_duplicates('key').set_index('key')

sel=pd.read_csv("coloc_egfr_overlap_loci.tsv",sep='\t')
print("loci to coloc:",len(sel),flush=True)
results=[]
for chrom in sorted(sel.chrom.unique()):
    pth=f"{G}/{chrom}/nominal/{chrom}.sig.cis_qtl_pairs.{chrom}.parquet"
    if chrom in ['chr1','chr2']: pth=f"{G}/{chrom}/nominal/{chrom}.sig.merged.parquet"
    if not os.path.exists(pth): continue
    sub=sel[sel.chrom==chrom]; chrn=chrom.replace('chr','')
    want=set(sub.best_cpg)
    pf=pq.ParquetFile(pth); frames=[]
    for i in range(pf.num_row_groups):
        t=pf.read_row_group(i,columns=['phenotype_id','variant_id','af','slope','slope_se']).to_pandas()
        frames.append(t[t.phenotype_id.isin(want)])
    md=pd.concat(frames) if frames else pd.DataFrame()
    for _,row in sub.iterrows():
        m=md[md.phenotype_id==row.best_cpg].copy()
        if len(m)<20: continue
        p=m.variant_id.str.split(':',expand=True)
        m['pos38']=p[1].astype(int); m['ref']=p[2].str.upper(); m['alt']=p[3].str.upper()
        cpos=int(row.best_cpg.split('_')[1])
        m=m[(m.pos38>cpos-100000)&(m.pos38<cpos+100000)]
        if len(m)<20: continue
        m['pos19']=[(lambda r:r[0][1]+1 if r else np.nan)(lo.convert_coordinate(chrom,x-1)) for x in m.pos38]
        m=m.dropna(subset=['pos19']); m['pos19']=m.pos19.astype(int)
        m['key']=chrn+':'+m.pos19.astype(str)
        mm=m[m.key.isin(gw.index)]
        if len(mm)<20: continue
        g=gw.loc[mm.key]
        sign=np.where(mm.alt.values==g.Allele1.str.upper().values,1,np.where(mm.ref.values==g.Allele1.str.upper().values,-1,0))
        ok=sign!=0
        if ok.sum()<20: continue
        b1=mm.slope.values[ok]; v1=(mm.slope_se.values[ok])**2
        b2=(g.Effect.values*sign)[ok]; v2=(g.StdErr.values[ok])**2
        good=np.isfinite(b1)&np.isfinite(v1)&np.isfinite(b2)&np.isfinite(v2)&(v1>0)&(v2>0)
        if good.sum()<20: continue
        post=coloc_abf(b1[good],v1[good],b2[good],v2[good])
        z1=b1[good]/np.sqrt(v1[good]); z2=b2[good]/np.sqrt(v2[good])
        zc=np.corrcoef(np.abs(z1),np.abs(z2))[0,1]
        results.append((row.best_cpg,row.nearest_gene,int(good.sum()),round(post[4],3),round(post[3],3),round(zc,3),round(row.asm_frac_fisher,3)))
res=pd.DataFrame(results,columns=['cpg','gene','nsnp','PP4','PP3','z_corr','asm']).sort_values('PP4',ascending=False)
res.to_csv("coloc_eGFR_results.tsv",sep='\t',index=False)
print("\n=== COLOCALIZING with eGFR (PP4>0.8) ===")
print(res[res.PP4>0.8].head(30).to_string(index=False))
print("\ntested:",len(res),"| PP4>0.8:",int((res.PP4>0.8).sum()),"| PP4>0.5:",int((res.PP4>0.5).sum()))
print("suggestive (z_corr>0.7 & PP3>0.5):",int(((res.z_corr>0.7)&(res.PP3>0.5)).sum()))
