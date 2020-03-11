include: "get_data.smk"
configfile: "conf/snakemake_config.json"


#RUNS, = glob_wildcards("mgf/{run}.mgf")
RUNS = get_runs(config["download"]["pxd_identifier"], ['raw'], config["download"]["file_pattern"])


rule search_targets:
	input:
		expand("mzid/{run}.pout", run=RUNS),
		expand("mzid/{run}.pout_dec", run=RUNS)


rule run_msgfplus:
	input:
		spectrum_file="mgf/{run}.mgf",
		msgfplus_conf=config["search"]["msgfplus_conf"],
		fasta=config["search"]["fasta"]
	output:
		"mzid/{run}.mzid"
	log:
		"logs/msgfplus/{run}.log"
	threads: config['search']['threads_per_search']
	shell:
		"""
		{config[search][msgfplus_exec]} -thread '{config[search][threads_per_search]}' -conf '{input.msgfplus_conf}' -d '{input.fasta}' -s '{input.spectrum_file}' -addFeatures 1
		mkdir -p mzid
		mv 'mgf/{wildcards.run}.mzid' 'mzid/{wildcards.run}.mzid'
		"""


rule create_pin:
	input:
		"mzid/{run}.mzid"
	output:
		"mzid/{run}.pin"
	log:
		"logs/msgf2pin/{run}.log"
	shell:
		"msgf2pin -P XXX '{input}' > '{output}'"


rule run_percolator:
	input:
		"mzid/{run}.pin"
	output:
		pout="mzid/{run}.pout",
		pout_dec="mzid/{run}.pout_dec"
	log:
		"logs/percolator/{run}.log"
	shell:
		"percolator --post-processing-tdc -U -m '{output.pout}' -M '{output.pout_dec}' '{input}'"
