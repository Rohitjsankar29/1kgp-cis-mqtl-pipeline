# NCI/Gadi run order

## 0. Clone and configure

```bash
cd /scratch/cy94/$USER
git clone https://github.com/RohitJayaSankar29/1kgp-cis-mqtl-pipeline.git
cd 1kgp-cis-mqtl-pipeline/1kgp_cis_mqtl_pipeline
nano config/gadi.env
```

## 1. Smoke test TensorQTL

```bash
python scripts/tests/00_test_tensorqtl.py --test-dir /scratch/cy94/$USER/test_tensorqtl
```

## 2. Build manifest

```bash
python scripts/manifest/00_build_manifest.py
```

or for a fixed sample list:

```bash
python scripts/manifest/00_make_manifest.py \
  --samples-tsv config/samples.tsv \
  --out config/samples_10.tsv \
  --n-samples 10
```

## 3. Extract methylation

Use the PBS array approach. Start very small.

```bash
qsub -J 1-1 scripts/methylation/02_stream_bam_extract_methylation.sh
```

Then scale to 2 samples, then 10.

## 4. Aggregate methylation matrix

```bash
python scripts/matrix/01_aggregate_bedmethyl.py \
  --input-dir /g/data/cy94/$USER/1kgp-mqtl/methylation/bedmethyl \
  --pattern '*_combined.bedmethyl.gz' \
  --output-dir /g/data/cy94/$USER/1kgp-mqtl/matrix \
  --min-coverage 5 \
  --min-samples 0.8
```

## 5. Prepare genotype PLINK files

```bash
bash scripts/genotype/04_prepare_plink_chr.sh chr22 input.vcf.gz /g/data/cy94/$USER/1kgp-mqtl/plink
```

## 6. Build covariates

```bash
python scripts/covariates/02_prepare_covariates.py \
  --samples config/samples.tsv \
  --methylation-pcs /g/data/cy94/$USER/1kgp-mqtl/matrix/methylation_beta.bed.gz \
  --geno-pcs /g/data/cy94/$USER/1kgp-mqtl/plink/chr22.eigenvec \
  --output /g/data/cy94/$USER/1kgp-mqtl/covariates/covariates.tsv
```

## 7. Run TensorQTL cis-mQTL

```bash
python scripts/tensorqtl/05_run_tensorqtl_pilot.py \
  --plink-prefix /g/data/cy94/$USER/1kgp-mqtl/plink/chr22 \
  --phenotype-bed /g/data/cy94/$USER/1kgp-mqtl/matrix/methylation_Mval.bed.gz \
  --covariates /g/data/cy94/$USER/1kgp-mqtl/covariates/covariates.tsv \
  --nominal-out /g/data/cy94/$USER/1kgp-mqtl/tensorqtl/nominal \
  --permutation-out /g/data/cy94/$USER/1kgp-mqtl/tensorqtl/permutation/chr22.permutation.txt.gz \
  --cis-window 1000000 \
  --n-permutations 1000 \
  --region chr22 \
  --skip-nominal
```

## 8. QC results

```bash
python scripts/qc/06_results_qc.py \
  --permutation /g/data/cy94/$USER/1kgp-mqtl/tensorqtl/permutation/chr22.permutation.txt.gz \
  --nominal-dir /g/data/cy94/$USER/1kgp-mqtl/tensorqtl/nominal \
  --out-dir /g/data/cy94/$USER/1kgp-mqtl/qc \
  --label pilot.chr22
```

## 9. Fine-mapping, SV, annotation, prioritisation

Run these only after TensorQTL outputs are verified.

```bash
Rscript scripts/finemapping/08_run_susie_finemap.R --help
python scripts/sv/09_integrate_sv_mqtl.py --help
python scripts/annotation/10_annotate_functional.py --help
python scripts/prioritisation/11_score_prioritise_variants.py --help
```
