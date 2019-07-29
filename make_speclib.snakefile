include: "get_data.snakefile"
configfile: "conf/snakemake_config.json"


RUNS, = glob_wildcards("mzid/{run}.pout")
#RUNS = get_runs(config["download"]["pxd_identifier"], ['raw'], config["download"]["file_pattern"])


rule speclib_targets:
    input:
        "speclib/spectral_library.peprec",
        "speclib/spectral_library.mgf"


rule run_pout_to_speclib:
    input:
        expand("mzid/{run}.pout", run=RUNS),
        expand("mgf/{run}.mgf", run=RUNS)
    output:
        "speclib/spectral_library.peprec",
        "speclib/spectral_library.mgf",
        temp(expand("mzid/{run}.pout_fixed", run=RUNS))
    log:
        "logs/pout_to_speclib/log.log"
    shell:
        """
        python3 scripts/pout_to_speclib.py -c conf/snakemake_config.json -i {config[download][pxd_identifier]} -p mzid -m mgf -o speclib -t 0.01
        """
