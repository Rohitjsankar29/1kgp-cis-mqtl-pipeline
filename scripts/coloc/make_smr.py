import sys, pandas as pd, numpy as np, os, pyarrow.parquet as pq
G="/g/data/cy94/rs4477/downstream/genome_452"
OUT="/scratch/cy94/rs4477/coloc/smr_input"; os.makedirs(OUT,exist_ok=True)
CHROM=sys.argv[1]
cands=set(l.strip() for l in open("/scratch/cy94/rs4477/coloc/smr_candidate_cpgs.txt") if l.strip().startswith(CHROM+"_"))
if not cands: print(f"{CHROM}: no candidate CpGs"); sys.exit(0)
pth=f"{G}/{CHROM}/nominal/{CHROM}.sig.cis_qtl_pairs.{CHROM}.parquet"
if CHROM in ['chr1','chr2']: pth=f"{G}/{CHROM}/nominal/{CHROM}.sig.merged.parquet"
chrn=CHROM.replace('chr','')
pf=pq.ParquetFile(pth); frames=[]
for i in range(pf.num_row_groups):
    d=pf.read_row_group(i,columns=['phenotype_id','variant_id','af','slope','slope_se','pval_nominal']).to_pandas()
    frames.append(d[d.phenotype_id.isin(cands)])
m=pd.concat(frames) if frames else pd.DataFrame()
print(f"{CHROM}: {m.phenotype_id.nunique()} CpGs, {len(m)} pairs",flush=True)
p=m.variant_id.str.split(':',expand=True)
m['pos']=p[1].astype(int); m['A1']=p[3].str.upper(); m['A2']=p[2].str.upper()  # alt=A1(effect), ref=A2
m['cpg_pos']=m.phenotype_id.str.split('_').str[1].astype(int)
# SNP ID = the variant_id (matches your 452 plink .bim exactly), all hg38
out=pd.DataFrame({'SNP':m.variant_id,'Chr':chrn,'BP':m.pos,'A1':m.A1,'A2':m.A2,'Freq':m.af.round(4),
    'Probe':m.phenotype_id,'Probe_Chr':chrn,'Probe_bp':m.cpg_pos,'Gene':m.phenotype_id,'Orientation':'+',
    'b':m.slope,'se':m.slope_se,'p':m.pval_nominal})
out.to_csv(f"{OUT}/{CHROM}.query.txt",sep='\t',index=False)
print(f"{CHROM}: wrote {len(out)} rows (hg38, no liftover)",flush=True)
