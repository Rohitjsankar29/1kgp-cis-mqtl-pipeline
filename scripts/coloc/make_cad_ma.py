import pandas as pd, numpy as np
from pyliftover import LiftOver
lo=LiftOver('hg19','hg38')
print("loading CAD primary (3GB, may take a bit)...",flush=True)
g=pd.read_csv('/scratch/cy94/rs4477/coloc/gwas/cad/CAD_GWAS_primary_discovery_meta.tsv',sep='\t',
              usecols=['CHR','BP','Allele1','Allele2','Freq1','Effect','StdErr','P-value','N'])
g=g.dropna(subset=['CHR','BP','Effect','StdErr','P-value'])
g['A1']=g.Allele1.str.upper(); g['A2']=g.Allele2.str.upper()
print("lifting b37->hg38...",flush=True)
def lift(c,p):
    r=lo.convert_coordinate('chr'+str(c),int(p)-1); return r[0][1]+1 if r else 0
g['pos38']=[lift(c,p) for c,p in zip(g.CHR,g.BP)]
g=g[g.pos38>0].copy()
g['SNP']='chr'+g.CHR.astype(str)+':'+g.pos38.astype(str)+':'+g.A2+':'+g.A1
g['SNP_rev']='chr'+g.CHR.astype(str)+':'+g.pos38.astype(str)+':'+g.A1+':'+g.A2
out=g.rename(columns={'Freq1':'freq','Effect':'b','StdErr':'se','P-value':'p'})[['SNP','SNP_rev','A1','A2','freq','b','se','p','N']]
out.to_csv('/scratch/cy94/rs4477/coloc/gwas/cad_hg38.ma.full',sep='\t',index=False)
print('CAD lifted:',len(out),flush=True)
