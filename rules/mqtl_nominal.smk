# =============================================================================
# mqtl_nominal  —  WORKED TEMPLATE (converted from run_nominal_gw.pbs)
# Copy this pattern for every other stage.
#
# What the .pbs did and where it went:
#   #PBS headers ............ deleted  -> become `resources:` (the PBS profile
#                                         turns them back into qsub flags)
#   micromamba activate ..... kept     -> runs inside the rule's shell
#   for N in seq 1 22 ....... deleted  -> the {chrom} wildcard + expand() below;
#                                         Snakemake runs the 22 chromosomes as
#                                         SEPARATE jobs, in parallel
#   [ -s $C.sig.done ] skip . deleted  -> Snakemake skips any chrom whose output
#                                         already exists (free resume)
#
# CANONICAL OUTPUT:
#   Big chromosomes (chr1, chr2) were run split-in-halves then merged, so their
#   final file is  {chrom}.sig.merged.parquet,  while chr3..22 produced a single
#   {chrom}.sig.cis_qtl_pairs.{chrom}.parquet. To give the pipeline ONE stable
#   name per chrom, the rule emits {chrom}.nominal.parquet as a symlink to
#   whichever real file exists. Downstream rules depend on {chrom}.nominal.parquet
#   and never have to know how it was computed.
# =============================================================================

B   = config["mqtl_in"]       # was  B=/scratch/cy94/rs4477/1kgp-cis-mqtl
G   = config["downstream"]    # was  G=/g/data/cy94/rs4477/downstream/genome_452
ENV = config["env"]
PM  = config["params"]
CHROMS = [f"chr{i}" for i in range(1, 23)]


rule mqtl_nominal:
    """cis-mQTL nominal pass for ONE chromosome (the body of the old for-loop)."""
    input:
        plink_bed = f"{B}/genotypes/genome_452/{{chrom}}/{{chrom}}.452samples.bed",
        pheno     = f"{B}/matrix/genome_452_cpgfilt/{{chrom}}.methylation_Mval.cpgfilt.bed.gz",
        cov       = f"{B}/covariates/genome_452/{{chrom}}.covariates.tsv",
        sigcpgs   = f"{B}/tensorqtl/genome_452/sig_by_chrom/{{chrom}}.sig_cpgs.txt",
    output:
        # one stable name per chrom; a symlink to the real parquet (see shell)
        nominal = f"{G}/{{chrom}}/nominal/{{chrom}}.nominal.parquet",
    params:
        plink_prefix = lambda w: f"{B}/genotypes/genome_452/{w.chrom}/{w.chrom}.452samples",
        outdir       = lambda w: f"{G}/{w.chrom}/nominal",
        prefix       = lambda w: f"{w.chrom}.sig",
        script       = "scripts/cov_v2/05c_run_nominal_significant.py",
        cis_window   = PM["cis_window"],
        maf          = PM["maf_threshold"],
        mamba        = ENV["mamba_exe"],
        mamba_root   = ENV["mamba_root"],
        env          = ENV["tensorqtl_env"],
    resources:
        queue    = "normal",     # your real queue was `normal` (CPU), NOT gpuvolta
        ncpus    = 12,
        mem      = "48GB",
        walltime = "08:00:00",
        gpu      = 0,
    threads: 12
    shell:
        r"""
        eval "$({params.mamba} shell hook --shell bash --root-prefix {params.mamba_root})"
        micromamba activate {params.env}
        mkdir -p {params.outdir}

        # only run the mapping if neither real output exists yet
        MERGED={params.outdir}/{wildcards.chrom}.sig.merged.parquet
        SINGLE={params.outdir}/{wildcards.chrom}.sig.cis_qtl_pairs.{wildcards.chrom}.parquet
        if [ ! -s "$MERGED" ] && [ ! -s "$SINGLE" ]; then
            python {params.script} \
              --plink-prefix {params.plink_prefix} \
              --phenotype-bed {input.pheno} \
              --covariates {input.cov} \
              --significant-cpgs {input.sigcpgs} \
              --out-dir {params.outdir} \
              --prefix {params.prefix} \
              --cis-window {params.cis_window} \
              --maf-threshold {params.maf} \
              --inverse-normal \
              --threads {threads}
        fi

        # point the canonical name at whichever real file exists (relative symlink)
        if   [ -s "$MERGED" ]; then ln -sf {wildcards.chrom}.sig.merged.parquet {output.nominal}
        elif [ -s "$SINGLE" ]; then ln -sf {wildcards.chrom}.sig.cis_qtl_pairs.{wildcards.chrom}.parquet {output.nominal}
        else echo "ERROR: no nominal parquet produced for {wildcards.chrom}" >&2; exit 1
        fi
        """


rule mqtl_nominal_all:
    """Aggregator: asking for this fans out into 22 independent mqtl_nominal jobs.
       Run just this stage with:  snakemake mqtl_nominal_all --profile profile/pbs"""
    input:
        expand(f"{G}/{{chrom}}/nominal/{{chrom}}.nominal.parquet", chrom=CHROMS),


# =============================================================================
# config/config.yaml block (already added):
#   mqtl_in / downstream / env{mamba_exe,mamba_root,tensorqtl_env} / params{cis_window,maf_threshold}
# Snakefile:   include: "rules/mqtl_nominal.smk"
# =============================================================================
