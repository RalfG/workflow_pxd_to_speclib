include: "get_data.snakefile"
include: "search_data.snakefile"


RUNS = get_runs(config["download"]["pxd_identifier"], ['raw'], config["download"]["file_pattern"])


rule targets:
	input:
		expand("mzid/{run}.pout", run=RUNS),
		expand("mzid/{run}.pout_dec", run=RUNS)
