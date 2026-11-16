# =============================================================================
# cis-mQTL convergence framework — top-level Snakefile
# Run:  snakemake -n                      (dry-run the whole DAG)
#       snakemake --profile profile/pbs   (submit each rule as a qsub job)
# =============================================================================
configfile: "config/config.yaml"

# ---- shared globals used by every rule file ----
B      = config["mqtl_in"]                       # /scratch/.../1kgp-cis-mqtl
G      = config["downstream"]                    # /g/data/.../genome_452
ENV    = config["env"]
PM     = config["params"]
CHROMS = [f"chr{i}" for i in range(1, 23)]

include: "rules/mqtl_nominal.smk"
include: "rules/downstream.smk"

# ---- default target: the framework's terminal product -----------------------
rule all:
    input:
        f"{G}/genome.prioritised_variant_level.tsv",
# Branch targets (run explicitly): prioritise_all, calibration_all, ancestry_all
