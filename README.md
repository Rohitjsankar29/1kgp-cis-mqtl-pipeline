# 1KGP cis-mQTL CpG-centric prioritisation pipeline

This repository contains a merged NCI/Gadi-ready pipeline for developing a reproducible CpG-centric cis-mQTL prioritisation framework.



## Core design

The main rule is: never store the cohort BAMs permanently. Each BAM is staged to `$PBS_JOBFS`, processed with Modkit, compressed bedMethyl outputs are moved to `/g/data`, and the BAM is deleted.

## Recommended workflow

```text
00_test_tensorqtl.py
        ↓
00_build_manifest.py / make_manifest_parallel.py
        ↓
02_stream_bam_extract_methylation.sh
        ↓
01_aggregate_bedmethyl.py
        ↓
04_prepare_plink_chr.sh
        ↓
02_prepare_covariates.py or 04_build_covariates.py
        ↓
05_run_tensorqtl_pilot.py
        ↓
06_results_qc.py
        ↓
08_run_susie_finemap.R
        ↓
09_integrate_sv_mqtl.py
        ↓
10_annotate_functional.py
        ↓
11_score_prioritise_variants.py
```

## Folder layout

```text
config/                 Gadi/NCI environment paths
jobs/                   PBS job scripts
scripts/manifest/       BAM/sample manifest construction
scripts/methylation/    Modkit methylation extraction
scripts/matrix/         bedMethyl aggregation and methylation matrices
scripts/covariates/     genotype/methylation/categorical covariates
scripts/genotype/       VCF to PLINK preparation
scripts/tensorqtl/      cis-mQTL mapping
scripts/qc/             QQ, Manhattan, lambda GC and hit table QC
scripts/finemapping/    SuSiE fine-mapping
scripts/sv/             structural variant integration
scripts/annotation/     functional annotation overlap
scripts/prioritisation/ CpG-centric scoring
scripts/tests/          synthetic TensorQTL smoke test
scripts/prototype/      earlier simplified prototype scripts retained for reference
scripts/deprecated/     deprecated helper scripts retained for reference
```

## First Gadi test

Do not start with the full cohort. Start with one chromosome and a very small number of samples.

```bash
git clone https://github.com/RohitJayaSankar29/1kgp-cis-mqtl-pipeline.git
cd 1kgp-cis-mqtl-pipeline/1kgp_cis_mqtl_pipeline
```

Edit:

```bash
nano config/gadi.env
```

Then test TensorQTL first:

```bash
python scripts/tests/00_test_tensorqtl.py --test-dir /scratch/cy94/$USER/test_tensorqtl
```

Then run a tiny real pilot:

```text
1 sample × chr22
2 samples × chr22
2 samples × chr21 + chr22
```

Only scale after each stage produces the expected output files.

## Important outputs

```text
config/bam_manifest.tsv
methylation/bedmethyl/*_combined.bedmethyl.gz
methylation/bedmethyl/*_combined.bedmethyl.gz.tbi
matrix/methylation_beta.bed.gz
matrix/methylation_Mval.bed.gz
covariates/covariates.tsv
tensorqtl/permutation/*.txt.gz
tensorqtl/nominal/*.parquet
qc/*.png / *.pdf / *.tsv
prioritisation/prioritised_variant_cpg_pairs.tsv
```

## Safety notes

- Do not commit BAM, VCF, BED, parquet or result files to Git.
- Keep raw BAMs temporary only.
- Use `$PBS_JOBFS` for staged BAMs.
- Use `/g/data` for final compressed bedMethyl and results.
- Keep logs for every PBS job.
