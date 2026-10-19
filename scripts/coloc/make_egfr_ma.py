import pandas as pd, numpy as np
from pyliftover import LiftOver
lo=LiftOver('hg19','hg38')
print("loading eGFR...", flush=True)
g=pd.read_csv('/scratch/cy94/rs4477/coloc/gwas/eGFR_stanzick2021.gz',sep='\t',compression='gzip',
              usecols=['chr','pos','Allele1','Allele2','Freq1','Effect','StdErr','P.value','n'])
print("lifting b37->hg38...", flush=True)
def lift(c,p):
    r=lo.convert_coordinate('chr'+str(c),int(p)-1)
    return r[0][1]+1 if r else 0
g['pos38']=[lift(c,p) for c,p in zip(g.chr,g.pos)]
g=g[g.pos38>0].copy()
g['A1']=g.Allele1.str.upper(); g['A2']=g.Allele2.str.upper()
# mQTL SNP = chr:pos:ref:alt. Build both allele orders so we match either.
g['snp_fwd']='chr'+g.chr.astype(str)+':'+g.pos38.astype(str)+':'+g.A2+':'+g.A1
g['snp_rev']='chr'+g.chr.astype(str)+':'+g.pos38.astype(str)+':'+g.A1+':'+g.A2
out=g[['snp_fwd','snp_rev','A1','A2','Freq1','Effect','StdErr','P.value','n']].copy()
out.columns=['SNP','SNP_rev','A1','A2','freq','b','se','p','N']
out.to_csv('/scratch/cy94/rs4477/coloc/gwas/eGFR_hg38.ma.full',sep='\t',index=False)
print('eGFR lifted:', len(out), flush=True)
