import pandas as pd, numpy as np, glob, subprocess, os
from pyliftover import LiftOver
import pyarrow.parquet as pq
lo=LiftOver('hg38','hg19')
G="/g/data/cy94/rs4477/downstream/genome_452"
GWAS="/scratch/cy94/rs4477/coloc/gwas/eGFR_stanzick2021.gz"
# load eGFR once (all), index by chr:pos
print("loading eGFR...")
gw=pd.read_csv(GWAS,sep='\t',compression='gzip',usecols=['chr','pos','Allele1','Allele2','Freq1','Effect','StdErr','n'])
gw['key']=gw.chr.astype(str)+':'+gw.pos.astype(str)
gwi=gw.set_index('key')
# candidate loci: top ASM-corroborated, non-imprinted (take best CpG per locus)
cand=pd.read_csv("/scratch/cy94/rs4477/coloc/coloc_candidate_loci.tsv",sep='\t')
cand=cand[cand.asm_frac_fisher.fillna(0)>0.25].copy()  # strong ASM only
# get best_cpg per candidate variant
enr=pd.read_csv(f"{G}/genome.prioritised_variant_level.enriched.txt.gz",sep='\t',usecols=['variant_id','best_cpg'])
cand=cand.merge(enr,on='variant_id',how='left')
# unique CpGs to test (one coloc per CpG)
cpgs=cand[['best_cpg','chrom','nearest_gene']].dropna().drop_duplicates('best_cpg')
print("CpG loci to scan:", len(cpgs))
results=[]
# process per chromosome (open parquet once)
for chrom in sorted(cpgs.chrom.unique()):
    pf_path=f"{G}/{chrom}/nominal/{chrom}.sig.cis_qtl_pairs.{chrom}.parquet"
    if chrom in ['chr1','chr2']:
        pf_path=f"{G}/{chrom}/nominal/{chrom}.sig.merged.parquet"
    if not os.path.exists(pf_path): continue
    sub=cpgs[cpgs.chrom==chrom]
    chrn=chrom.replace('chr','')
    pf=pq.ParquetFile(pf_path)
    # load all needed CpGs for this chrom
    want=set(sub.best_cpg)
    frames=[]
    for i in range(pf.num_row_groups):
        t=pf.read_row_group(i,columns=['phenotype_id','variant_id','af','slope','slope_se']).to_pandas()
        frames.append(t[t.phenotype_id.isin(want)])
    md=pd.concat(frames) if frames else pd.DataFrame()
    for _,row in sub.iterrows():
        m=md[md.phenotype_id==row.best_cpg].copy()
        if len(m)<20: continue
        p=m.variant_id.str.split(':',expand=True)
        m['pos38']=p[1].astype(int); m['ref']=p[2]; m['alt']=p[3]
        # tight +/-100kb
        cpos=int(row.best_cpg.split('_')[1])
        m=m[(m.pos38>cpos-100000)&(m.pos38<cpos+100000)]
        if len(m)<20: continue
        m['pos19']=[ (lambda r: r[0][1]+1 if r else np.nan)(lo.convert_coordinate(chrom,x-1)) for x in m.pos38]
        m=m.dropna(subset=['pos19']); m['pos19']=m.pos19.astype(int)
        m['key']=chrn+':'+m.pos19.astype(str)
        m['maf']=np.where(m.af>0.5,1-m.af,m.af)
        mm=m[m.key.isin(gwi.index)]
        if len(mm)<20: continue
        g=gwi.loc[mm.key]
        # write coloc input
        sign=np.where(mm.alt.values==g.Allele1.str.upper().values,1,np.where(mm.ref.values==g.Allele1.str.upper().values,-1,0))
        ok=sign!=0
        df=pd.DataFrame({'snp':mm.key.values[ok],'mqtl_beta':mm.slope.values[ok],'mqtl_se':mm.slope_se.values[ok],'mqtl_maf':mm.maf.values[ok],'mqtl_n':452,'gwas_beta':(g.Effect.values*sign)[ok],'gwas_se':g.StdErr.values[ok],'gwas_maf':np.where(g.Freq1.values>0.5,1-g.Freq1.values,g.Freq1.values)[ok],'gwas_n':g.n.values[ok]}).dropna()
        if len(df)<20: continue
        tmp=f"/scratch/cy94/rs4477/coloc/tmp_scan"
        df.to_csv(tmp+".coloc_input.tsv",sep='\t',index=False)
        r=subprocess.run(['/scratch/cy94/rs4477/micromamba/envs/mqtl_tensorqtl/bin/Rscript','/scratch/cy94/rs4477/coloc/run_coloc.R',tmp],capture_output=True,text=True)
        for line in r.stdout.split('\n'):
            if 'PP4 (shared' in line:
                pp4=float(line.split(':')[1].strip())
                results.append((row.best_cpg,row.nearest_gene,len(df),pp4))
        print(f"{row.nearest_gene} {row.best_cpg}: done", flush=True)
res=pd.DataFrame(results,columns=['cpg','gene','nsnp','PP4'])
res=res.sort_values('PP4',ascending=False)
res.to_csv("/scratch/cy94/rs4477/coloc/coloc_eGFR_scan.tsv",sep='\t',index=False)
print("\n=== COLOCALIZING loci (PP4>0.8) with eGFR ===")
print(res[res.PP4>0.8].to_string(index=False))
print("\ntotal tested:",len(res),"| PP4>0.8:",int((res.PP4>0.8).sum()),"| PP4>0.5:",int((res.PP4>0.5).sum()))
