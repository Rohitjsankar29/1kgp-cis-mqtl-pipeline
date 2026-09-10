import sys, pandas as pd, numpy as np
from pyliftover import LiftOver
lo=LiftOver('hg19','hg38')
which=sys.argv[1]
def liftvec(chrs,poss):
    uniq={}; out=np.zeros(len(poss),dtype=np.int64)
    for i,(c,p) in enumerate(zip(chrs,poss)):
        k=(c,p)
        if k not in uniq:
            r=lo.convert_coordinate('chr'+str(c),int(p)-1); uniq[k]=r[0][1]+1 if r else 0
        out[i]=uniq[k]
    return out
if which=='t2d':
    print("T2D DIAMANTE-TA...",flush=True)
    g=pd.read_csv('/scratch/cy94/rs4477/coloc/gwas/diamante/DIAMANTE-EUR.sumstat.txt.gz',sep=' ',compression='gzip')
    g.columns=[c.strip() for c in g.columns]
    g=g.rename(columns={'chromosome(b37)':'chr','position(b37)':'pos','effect_allele':'A1','other_allele':'A2',
        'effect_allele_frequency':'freq','Fixed-effects_beta':'b','Fixed-effects_SE':'se',
        'Fixed-effects_p-value':'p'})
    g['N']=933970
    out='/scratch/cy94/rs4477/coloc/gwas/t2d_hg38.ma.full'
else:
    print("CAD primary...",flush=True)
    g=pd.read_csv('/scratch/cy94/rs4477/coloc/gwas/cad/CAD_GWAS_primary_discovery_meta.tsv',sep='\t',
        usecols=['CHR','BP','Allele1','Allele2','Freq1','Effect','StdErr','P-value','N'])
    g=g.rename(columns={'CHR':'chr','BP':'pos','Allele1':'A1','Allele2':'A2','Freq1':'freq',
        'Effect':'b','StdErr':'se','P-value':'p'})
    out='/scratch/cy94/rs4477/coloc/gwas/cad_hg38.ma.full'
g=g.dropna(subset=['chr','pos','b','se','p'])
g['A1']=g.A1.astype(str).str.upper(); g['A2']=g.A2.astype(str).str.upper()
print("lifting",len(g),"variants...",flush=True)
g['pos38']=liftvec(g.chr.values,g.pos.values)
g=g[g.pos38>0].copy()
g['SNP']='chr'+g.chr.astype(str)+':'+g.pos38.astype(str)+':'+g.A2+':'+g.A1
g['SNP_rev']='chr'+g.chr.astype(str)+':'+g.pos38.astype(str)+':'+g.A1+':'+g.A2
g[['SNP','SNP_rev','A1','A2','freq','b','se','p','N']].to_csv(out,sep='\t',index=False)
print("wrote",out,len(g),flush=True)
