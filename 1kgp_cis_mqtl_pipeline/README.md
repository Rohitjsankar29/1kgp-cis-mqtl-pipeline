# 1KGP cis-mQTL prioritisation pilot pipeline

This repository builds a CpG-centric cis-mQTL framework.

Workflow:
1. Build 1KGP ONT BAM manifest
2. Extract one chromosome at a time from modBAMs
3. Run modkit CpG 5mC methylation extraction
4. Aggregate bedMethyl into tensorQTL phenotype BED
5. Build covariates
6. Prepare genotype PCA / PLINK inputs
7. Run TensorQTL cis-mQTL
8. QC results
9. Build variant-CpG feature matrix for prioritisation

Storage rule: do not keep cohort BAMs. Stage or stream one sample/chromosome, keep only compressed bedMethyl and downstream matrices.
