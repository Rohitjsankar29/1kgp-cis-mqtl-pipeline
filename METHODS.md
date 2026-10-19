# Genome-wide cis-mQTL Prioritisation Framework - Methods

Reproducible CpG-centric cis-mQTL variant-prioritisation on 452 1000G-ONT samples.

## Cohort
452 unrelated founders (500 ONT roster intersect 2594 unrelated 1000G founders).
Sample list: config/keep_geno_452.txt

## Stages
1. Methylation: ONT modBAM -> modkit -> M-value matrix -> CpG-SNP filter (in-sample bim, ~2.7% dropped)
2. Genotypes: 1000G panel VCFs -> plink --snps-only, 452 samples (chr1/chr7 re-downloaded, bgzip corruption)
3. Covariates (24): sex + platform + 2 basecaller + 10 genoPC + 10 methPC
4. Discovery: tensorQTL permutation, 100 perms, 184 chunks -> 25.9M tested, 1.64M sig (FDR<0.05 BH)
5. Nominal: 05c_run_nominal_significant.py per-chromosome, mem ~1GB/1000 CpGs -> 205GB parquet
6. Fine-mapping: SuSiE (finemapsusie.py, tensorqtl.susie.map_susie summary mode, L=10) -> 4.25M cs variants, 192k PIP>0.9
7. SV: filter ALT symbolic OR SVTYPE; sv_annotate p<1e-3 AND r2>0.8 -> 29.1% CpGs implicated
8. Functional: functional_annotate.py per-chrom (nearest gene, CpG island, cCRE) -> 30.5% cCRE
9. Scoring: prioritise_framework.py (score_weighted + elasticnet) -> collapse -> enrich -> 4.25M ranked
10. ASM: read-level haplotype MWU (asm_readlevel_v2.py) - corroboration

Environment: environment/mqtl_tensorqtl.yml (tensorqtl 1.0.10, torch 2.4.1, susieR, rpy2)
Results: results/genome_452/headline_numbers.txt
Params: cis +/-1Mb, MAF 0.05, FDR BH, max L=10
