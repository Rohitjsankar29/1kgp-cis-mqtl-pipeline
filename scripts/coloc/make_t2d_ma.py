import pandas as pd, numpy as np
from pyliftover import LiftOver
lo=LiftOver('hg19','hg38')
print("loading DIAMANTE-TA...",flush=True)
g=pd.read_csv('/scratch/cy94/rs4477/coloc/gwas/diamante/DIAMANTE-TA.sumstat.txt.gz',sep=' ',compression='gzip')
g.columns=[c.strip() for c in g.columns]
g=g.rename(columns={'chromosome(b37)':'chr','position(b37)':'pos','effect_allele':'A1','other_allele':'A2',
                    'effect_allele_frequency':'freq','Fixed-effects_beta':'b','Fixed-effects_SE':'se',
                    'Fixed-effects_p-value':'p','Effective_sample_size':'N'})
g=g.dropna(subset=['chr','pos','b','se','p'])
g['A1']=g.A1.str.upper(); g['A2']=g.A2.str.upper()
print("lifting b37->hg38...",flush=True)
def lift(c,p):
    r=lo.convert_coordinate('chr'+str(c),int(p)-1); return r[0][1]+1 if r else 0
g['pos38']=[lift(c,p) for c,p in zip(g.chr,g.pos)]
g=g[g.pos38>0].copy()
g['SNP']='chr'+g.chr.astype(str)+':'+g.pos38.astype(str)+':'+g.A2+':'+g.A1
g['SNP_rev']='chr'+g.chr.astype(str)+':'+g.pos38.astype(str)+':'+g.A1+':'+g.A2
out=g[['SNP','SNP_rev','A1','A2','freq','b','se','p','N']]
out.to_csv('/scratch/cy94/rs4477/coloc/gwas/t2d_hg38.ma.full',sep='\t',index=False)
print('T2D lifted:',len(out),flush=True)
