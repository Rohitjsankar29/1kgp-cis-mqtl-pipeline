# Merge notes

This repository has been merged so that guide/Cursor-agent scripts are used for actual NCI execution, while the dissertation-specific scripts remain for downstream prioritisation.

## What changed

| Area | Merged choice |
|---|---|
| Manifest | Added guide `00_build_manifest.py`, `00_make_manifest.py`, `make_manifest_parallel.py`, `build_manifest_chr.py` |
| Methylation extraction | Added production `02_stream_bam_extract_methylation.sh` and single-sample `run_modkit_pileup.sh` |
| Matrix | Replaced simplified matrix script with guide `01_aggregate_bedmethyl.py`; retained pilot `03_aggregate_pilot.py` |
| Covariates | Added both general `02_prepare_covariates.py` and pilot `04_build_covariates.py` |
| TensorQTL | Replaced simplified script with guide `05_run_tensorqtl_pilot.py` |
| QC | Replaced simplified QC with guide `06_results_qc.py` |
| Testing | Added synthetic `00_test_tensorqtl.py` |
| Downstream | Kept SuSiE, SV, annotation and prioritisation scripts |
| Old scripts | Moved earlier simplified scripts into `scripts/prototype/` |

## Recommended execution scripts

Use these for actual NCI runs:

```text
scripts/manifest/00_build_manifest.py
scripts/methylation/02_stream_bam_extract_methylation.sh
scripts/matrix/01_aggregate_bedmethyl.py
scripts/covariates/02_prepare_covariates.py
scripts/tensorqtl/05_run_tensorqtl_pilot.py
scripts/qc/06_results_qc.py
scripts/finemapping/08_run_susie_finemap.R
scripts/sv/09_integrate_sv_mqtl.py
scripts/annotation/10_annotate_functional.py
scripts/prioritisation/11_score_prioritise_variants.py
```

## Scripts retained only for reference

```text
scripts/prototype/*
scripts/deprecated/geno_pca.py
```
