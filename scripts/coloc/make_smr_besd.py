# Convert mQTL nominal parquets -> SMR text format (for --make-besd)
# SMR "query" format columns: SNP Chr BP A1 A2 Freq Probe Probe_Chr Probe_bp Gene Orientation b se p
import sys, pandas as pd, numpy as np, glob, os, pyarrow.parquet as pq
from pyliftover import LiftOver
lo=LiftOver('hg38','hg19')  # SMR/GWAS in b37, mQTL in hg38
G="/g/data/cy94/rs4477/downstream/genome_452"
OUT="/scratch/cy94/rs4477/coloc/smr_input"
os.makedirs(OUT, exist_ok=True)
CHROM=sys.argv[1]
pth=f"{G}/{CHROM}/nominal/{CHROM}.sig.cis_qtl_pairs.{CHROM}.parquet"
if CHROM in ['chr1','chr2']: pth=f"{G}/{CHROM}/nominal/{CHROM}.sig.merged.parquet"
chrn=CHROM.replace('chr','')
pf=pq.ParquetFile(pth)
rows=[]
for i in range(pf.num_row_groups):
    d=pf.read_row_group(i,columns=['phenotype_id','variant_id','af','slope','slope_se','pval_nominal']).to_pandas()
    p=d.variant_id.str.split(':',expand=True)
    d['pos38']=p[1].astype(int); d['A1']=p[3]; d['A2']=p[2]  # alt=effect(A1), ref=A2
    # liftover to b37
    d['pos19']=[ (lambda r:r[0][1]+1 if r else 0)(lo.convert_coordinate(CHROM,x-1)) for x in d.pos38]
    d=d[d.pos19>0]
    # CpG probe position (hg38 -> b37 too)
    d['cpg_pos38']=d.phenotype_id.str.split('_').str[1].astype(int)
    d['cpg_pos19']=[ (lambda r:r[0][1]+1 if r else 0)(lo.convert_coordinate(CHROM,x-1)) for x in d.cpg_pos38]
    for _,r in d.iterrows():
        rows.append((r.variant_id, chrn, r.pos19, r.A1, r.A2, round(r.af,4),
                     r.phenotype_id, chrn, int(r.cpg_pos19), r.phenotype_id, '+',
                     r.slope, r.slope_se, r.pval_nominal))
o=pd.DataFrame(rows, columns=['SNP','Chr','BP','A1','A2','Freq','Probe','Probe_Chr','Probe_bp','Gene','Orientation','b','se','p'])
o.to_csv(f"{OUT}/{CHROM}.smr_query.txt", sep='\t', index=False)
print(f"{CHROM}: {len(o)} SNP-CpG pairs -> {OUT}/{CHROM}.smr_query.txt")
