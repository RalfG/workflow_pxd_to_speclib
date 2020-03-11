include: "get_data.smk"
configfile: "conf/snakemake_config.json"


#RUNS, = glob_wildcards("mzid/{run}.pout")
RUNS = get_runs(config["download"]["pxd_identifier"], ['raw'], config["download"]["file_pattern"])


rule speclib_targets:
    input:
        "speclib/calibrated_retention_times.peprec"


rule retention_time_calibration:
    input:
        expand("mzid/{run}.pout", run=RUNS),
        expand("mgf/{run}.mgf", run=RUNS)
    output:
        "speclib/calibrated_retention_times.peprec"
    conda:
        "envs/retention_time_calibration.yml"
    shell:
        """
        python scripts/retention_time_calibration.py --mgf="mgf" --mzid="mzid" --output-file="speclib/calibrated_retention_times.peprec"
        """
