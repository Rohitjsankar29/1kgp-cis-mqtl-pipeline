# chr22 5-sample cis-mQTL pilot results

## Samples
- HG02470
- HG02479
- HG03027
- HG03045
- HG03079

## Completed pilot stages

### 1. Methylation extraction
Per-sample chr22 bedMethyl files were generated using modkit.

Output directory:
/g/data/xl04/rs4477/1kgp-cis-mqtl/methylation/chr22/

### 2. Methylation matrix
5-sample chr22 methylation matrices were generated.

Outputs:
/g/data/xl04/rs4477/1kgp-cis-mqtl/matrix/chr22_5samples/methylation_beta.bed.gz
/g/data/xl04/rs4477/1kgp-cis-mqtl/matrix/chr22_5samples/methylation_Mval.bed.gz
/g/data/xl04/rs4477/1kgp-cis-mqtl/matrix/chr22_5samples/site_qc.tsv

TensorQTL-compatible phenotype:
/g/data/xl04/rs4477/1kgp-cis-mqtl/matrix/chr22_5samples_tensorqtl/methylation_Mval.bed.gz

### 3. Genotype input
chr22 1000 Genomes genotype VCF was downloaded and subset to the 5 methylation samples.

PLINK outputs:
/g/data/xl04/rs4477/1kgp-cis-mqtl/genotypes/chr22_5samples/chr22.5samples.bed
/g/data/xl04/rs4477/1kgp-cis-mqtl/genotypes/chr22_5samples/chr22.5samples.bim
/g/data/xl04/rs4477/1kgp-cis-mqtl/genotypes/chr22_5samples/chr22.5samples.fam

### 4. Covariates
Covariate file:
/g/data/xl04/rs4477/1kgp-cis-mqtl/covariates/chr22_5samples/covariates.tsv

### 5. TensorQTL technical test
TensorQTL successfully loaded the genotype, phenotype and covariate files. A technical test completed successfully.

Output directory:
/g/data/xl04/rs4477/1kgp-cis-mqtl/tensorqtl/chr22_5samples/tiny_test/

## Notes
Full nominal TensorQTL for chr22 was too resource-heavy for this 5-sample pilot. Several nominal chunk attempts produced memory, pyarrow, or duplicate-output issues. For this pilot, the successful TensorQTL technical test should be treated as framework validation, not biological inference.
