# =============================================================================
# downstream.smk  —  everything after nominal mapping.
# Dependency spine:
#   sv_extract -> sv_aggregate ┐
#   split_beds -> functional   ├-> prioritise -> combine -> collapse  (THE FRAMEWORK)
#   finemap (EXTERNAL) ────────┘
# Parallel branches: coloc (disease), calibration (QC), confirm (1000-perm),
#                    ancestry (SMR sensitivity), asm (per-sample).
#
# TAGS you'll see:
#   # EXTERNAL = input produced by a stage whose script isn't in jobs/genome_452
#               (finemap SuSiE, the 100-perm sig-CpG discovery, ASM rescore).
#               Treated as pre-existing; wire its own rule when you locate it.
#   # VERIFY   = output filename inferred from the .pbs; confirm with `ls` and
#               fix if wrong (same one-line check we did for nominal).
# =============================================================================

REPO = config["repo"]                       # /scratch/.../1kgp-cis-mqtl-pipeline
SC   = config["scripts_cov2"]               # /scratch/cy94/rs4477/cov_v2
A    = config["annotation"]                 # /scratch/.../reference/annotation
COLOC = config["coloc_dir"]                 # /scratch/cy94/rs4477/coloc
EROOT = f'{ENV["mamba_root"]}/envs'
PY    = f'{EROOT}/{ENV["tensorqtl_env"]}/bin/python3'
BCF   = f'{EROOT}/mqtl_modkit/bin/bcftools'
TABIX = f'{EROOT}/mqtl_modkit/bin/tabix'
SUPERPOPS = ["EUR", "AFR", "SAS", "EAS", "AMR"]


# ---- shared: per-chromosome annotation BED splits (was awk inside two jobs) --
rule split_beds:
    input:
        tss    = f"{A}/gencode_tss_genome.bed",
        island = f"{A}/cpg_islands_genome.bed",
        ccre   = f"{A}/ccre_4col.bed",
    output:
        tss    = f"{A}/perchrom/{{chrom}}.tss.bed",
        island = f"{A}/perchrom/{{chrom}}.island.bed",
        ccre   = f"{A}/perchrom/{{chrom}}.ccre.bed",
    shell:
        r"""
        awk -v c={wildcards.chrom} '$1==c' {input.tss}    > {output.tss}
        awk -v c={wildcards.chrom} '$1==c' {input.island} > {output.island}
        awk -v c={wildcards.chrom} '$1==c' {input.ccre}   > {output.ccre}
        """


# ---- SV axis: extract SV genotypes, then aggregate to per-variant linkage ----
rule sv_extract:                                   # <- sv_extract_gw.pbs
    input:
        panel = f"{B}/genotypes/{{chrom}}/1kGP_high_coverage_Illumina.{{chrom}}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz",
        keep  = f"{B}/config/keep_geno_452.txt",
    output:
        gt = f"{B}/sv_genome_452/{{chrom}}.sv.gt.tsv",
    resources: queue="normal", ncpus=4, mem="16GB", walltime="05:00:00", gpu=0
    shell:
        r"""
        {BCF} view -r {wildcards.chrom} -S {input.keep} --force-samples {input.panel} \
          | {BCF} view -i 'ALT[0]~"<" || INFO/SVTYPE!="."' -Oz -o {output.gt}.vcf.gz
        {BCF} index -t {output.gt}.vcf.gz
        {{
          printf "sv_id\tchrom\tpos\tend\tsvtype\tsvlen"
          for s in $({BCF} query -l {output.gt}.vcf.gz); do printf "\t%s" "$s"; done
          printf "\n"
          {BCF} query -f '%ID\t%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%INFO/SVLEN[\t%GT]\n' {output.gt}.vcf.gz
        }} > {output.gt}
        """

rule sv_aggregate:                                 # <- sv_agg.pbs / sv_aggregate.py
    # argless script: reads all per-chrom sv.gt.tsv, writes per-chrom sv_annotation
    input:
        expand(f"{B}/sv_genome_452/{{chrom}}.sv.gt.tsv", chrom=CHROMS),
    output:
        # VERIFY: prioritise consumes $G/$C/sv/$C.sv_annotation.txt.gz
        expand(f"{G}/{{chrom}}/sv/{{chrom}}.sv_annotation.txt.gz", chrom=CHROMS),
    resources: queue="normal", ncpus=12, mem="64GB", walltime="02:00:00", gpu=0
    shell:
        r"{PY} {SC}/sv_aggregate.py"


# ---- functional axis: TSS / CpG-island / cCRE annotation per chrom ----------
rule functional:                                   # <- func_annot_gw.pbs
    input:
        matrix  = f"{B}/matrix/genome_452_cpgfilt/{{chrom}}.methylation_Mval.cpgfilt.bed.gz",
        sigcpgs = f"{B}/tensorqtl/genome_452/sig_by_chrom/{{chrom}}.sig_cpgs.txt",   # EXTERNAL (100-perm discovery)
        tss     = rules.split_beds.output.tss,
        island  = rules.split_beds.output.island,
        ccre    = rules.split_beds.output.ccre,
    output:
        f"{G}/{{chrom}}/functional/{{chrom}}.functional.tsv",
    resources: queue="normal", ncpus=8, mem="32GB", walltime="03:00:00", gpu=0
    shell:
        r"""
        {PY} {SC}/functional_annotate.py \
          --matrix {input.matrix} --significant-cpgs {input.sigcpgs} \
          --tss-bed {input.tss} --cpg-island-bed {input.island} --ccre-bed {input.ccre} \
          --out {output}
        """


# ---- prioritise: the 5-axis convergence score (UNIFIES gw + c12 + big) ------
rule prioritise:                                   # <- prioritise_gw/_c12/_big.pbs (now ONE rule)
    input:
        susie      = f"{G}/{{chrom}}/finemap/{{chrom}}.susie.tsv",            # EXTERNAL (SuSiE fine-map)
        bim        = f"{B}/genotypes/genome_452/{{chrom}}/{{chrom}}.452samples.bim",
        sv         = f"{G}/{{chrom}}/sv/{{chrom}}.sv_annotation.txt.gz",
        functional = rules.functional.output,
        nominal    = f"{G}/{{chrom}}/nominal/{{chrom}}.nominal.parquet",       # canonical (from mqtl_nominal.smk)
        tss        = rules.split_beds.output.tss,
        ccre       = rules.split_beds.output.ccre,
        island     = rules.split_beds.output.island,
    output:
        f"{G}/{{chrom}}/{{chrom}}.prioritised.txt.gz",
    resources: queue="normal", ncpus=8, mem="48GB", walltime="05:00:00", gpu=0
    shell:
        r"""
        {PY} {SC}/prioritise_framework.py \
          --susie {input.susie} --bim {input.bim} --sv {input.sv} \
          --functional {input.functional} --nominal {input.nominal} \
          --tss-bed {input.tss} --ccre-bed {input.ccre} --cpg-island-bed {input.island} \
          --out {output}
        """


# ---- combine per-chrom -> genome table, then collapse to variant level -------
rule combine_prioritised:                          # (concat step; was done by hand)
    input:
        expand(f"{G}/{{chrom}}/{{chrom}}.prioritised.txt.gz", chrom=CHROMS),
    output:
        f"{G}/genome.prioritised.tsv",
    run:
        import pandas as pd
        frames = [pd.read_csv(f, sep="\t", compression="infer") for f in input]
        pd.concat(frames, ignore_index=True).to_csv(output[0], sep="\t", index=False)

rule collapse:                                     # <- collapse.pbs
    input:
        rules.combine_prioritised.output,
    output:
        f"{G}/genome.prioritised_variant_level.tsv",
    resources: queue="normal", ncpus=12, mem="64GB", walltime="01:00:00", gpu=0
    shell:
        r"""
        {PY} {SC}/collapse_to_variant_level.py \
          --prioritised {input} --out {output} --score score_weighted
        """
# NOTE: your thesis's final table is genome.prioritised_final_asm_rescored.txt.gz,
# i.e. there's an ASM-rescore step (rescore_with_asm.py) between prioritise and
# final that isn't in jobs/genome_452. Wire it as a rule between combine and
# collapse once you locate that script.  # EXTERNAL


# ===== PARALLEL BRANCHES (not on the prioritisation spine) ====================

# ---- disease integration: coloc per chrom (eGFR) ----------------------------
rule coloc:                                        # <- coloc_egfr_array.pbs / fast_coloc_v2.py
    output:
        done = touch(f"{COLOC}/chr{{chrom_n}}.coloc.done"),   # VERIFY: replace with real coloc output file
    resources: queue="normal", ncpus=2, mem="24GB", walltime="05:00:00", gpu=0
    shell:
        r"cd {COLOC} && {PY} fast_coloc_v2.py chr{wildcards.chrom_n}"

# ---- calibration: genome-wide null permutation (GPU) + lambda ---------------
rule calibration_null:                             # <- lambda_null_genomewide.pbs (GPU)
    input:
        plink = f"{B}/genotypes/genome_452/{{chrom}}/{{chrom}}.452samples.bed",
        pheno = f"{B}/matrix/genome_452_cpgfilt/{{chrom}}.methylation_Mval.cpgfilt.bed.gz",
        cov   = f"{B}/covariates/genome_452/{{chrom}}.covariates.tsv",
    output:
        f"{G}/calibration/null_gw/{{chrom}}.null.perm.txt.gz",
    resources: queue="gpuvolta", ncpus=12, mem="90GB", walltime="02:00:00", gpu=1
    shell:
        r"""
        {PY} {SC}/05d_permuted_null.py \
          --plink-prefix {B}/genotypes/genome_452/{wildcards.chrom}/{wildcards.chrom}.452samples \
          --phenotype-bed {input.pheno} --covariates {input.cov} \
          --out {output} --cis-window 1000000 --maf-threshold 0.05 \
          --n-permutations 100 --n-phenotypes 10000 --inverse-normal --threads 12
        """

# ---- 1000-perm confirmation on prioritised hits -----------------------------
rule confirm_1000:                                 # <- confirm_perm_1000.pbs
    input:
        cpgs  = f"{REPO}/config/confirmation_cpgs.txt",
        pheno = f"{B}/matrix/genome_452_cpgfilt/{{chrom}}.methylation_Mval.cpgfilt.bed.gz",
        cov   = f"{B}/covariates/genome_452/{{chrom}}.covariates.tsv",
    output:
        f"{G}/{{chrom}}/confirm/{{chrom}}.confirm1000.txt.gz",
    resources: queue="normal", ncpus=12, mem="48GB", walltime="06:00:00", gpu=0
    shell:
        r"""
        LIST=$TMPDIR/{wildcards.chrom}.cpgs
        grep "^{wildcards.chrom}_" {input.cpgs} > $LIST || true
        BED=$TMPDIR/{wildcards.chrom}.subset.bed
        zcat {input.pheno} | awk -v L=$LIST 'BEGIN{{while((getline x<L)>0)k[x]=1}} NR==1||($4 in k)' > $BED
        gzip -f $BED
        {PY} {SC}/05a_run_permutation.py \
          --plink-prefix {B}/genotypes/genome_452/{wildcards.chrom}/{wildcards.chrom}.452samples \
          --phenotype-bed $BED.gz --covariates {input.cov} \
          --out {output} --n-permutations 1000 --inverse-normal --threads 12
        """

# ---- ancestry SMR sensitivity, per superpopulation --------------------------
rule ancestry_sensitivity:                         # <- ancestry_eur.pbs / ancestry_persuperpop.pbs
    output:
        f"{G}/ancestry/sensitivity/ancestry_sensitivity_summary_{{pop}}.tsv",
    resources: queue="normal", ncpus=8, mem="32GB", walltime="03:00:00", gpu=0
    shell:
        r"""
        bash {REPO}/scripts/ancestry/02_compute_ancestry_ld.sh {wildcards.pop}
        bash {REPO}/scripts/ancestry/03_run_smr_ancestry.sh    {wildcards.pop}
        {PY} {REPO}/scripts/ancestry/04_compare_ancestry_vs_multi.py \
          --pop {wildcards.pop} --primary-dir {COLOC} \
          --ancestry-dir {G}/ancestry/smr/{wildcards.pop} \
          --out-dir {G}/ancestry/sensitivity
        """


# ===== aggregators (ask for these as targets) =================================
rule prioritise_all:
    input: expand(f"{G}/{{chrom}}/{{chrom}}.prioritised.txt.gz", chrom=CHROMS)

rule calibration_all:
    input: expand(f"{G}/calibration/null_gw/{{chrom}}.null.perm.txt.gz", chrom=CHROMS)

rule ancestry_all:
    input: expand(f"{G}/ancestry/sensitivity/ancestry_sensitivity_summary_{{pop}}.tsv", pop=SUPERPOPS)
