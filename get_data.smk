configfile: "conf/snakemake_config.json"


from scripts.download_pride_project import get_files_df

def get_runs(pxd_identifier, extensions, file_pattern):
	runs =  list(get_files_df(pxd_identifier,
		extensions,
		file_pattern
	)['fileName'].str.replace('.raw', '', case=False))
	return runs

RUNS = get_runs(config["download"]["pxd_identifier"], ['raw'], config["download"]["file_pattern"])


rule download_targets:
	input:
		expand("mgf/{run}.mgf", run=RUNS)


rule download:
	input:
	output:
		expand("raw/{run}.raw", run=RUNS)
	log:
		"logs/download_pride_project/log.log"
	shell:
		"python3 scripts/download_pride_project.py -f raw -p '{config[download][file_pattern]}' '{config[download][pxd_identifier]}'"


rule convert_to_mgf:
	input:
		"raw/{run}.raw"
	output:
		"mgf/{run}.mgf"
	shell:
		"{config[convert][exec]} --input='{input}' --output_file='{output}' -f=0 -m=0"
