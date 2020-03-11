# Workflow: pxd_to_speclib
Snakemake workflow that reanalyzes proteomics data from a PRIDE Archive project to create a spectral library.

The workflow goes through the following steps:
- Download RAW files from PRIDE Archive for given PXD identifier
- Convert RAW files to MGF using the [CompOmics ThermoRawFileParser](https://github.com/compomics/ThermoRawFileParser)
- Search with [MSGFPlus](https://github.com/MSGFPlus/msgfplus)
- Generate Percolator input files
- Postprocess search results with [Percolator](https://github.com/percolator/percolator/)
- Parse results to generate a spectral library

## Requirements
- Conda (tested on Linux)

## Run the workflow
1. Create and activate the environment:
```
conda env create -f envs/pxd_to_speclib.yml
conda activate pxd_to_speclib
```
2. Setup your configuration: 
    - `conf/snakemake_config.json` (see [Configuration](#configuration))
    - `conf/msgfplus_params.txt`
    - Add required input files (e.g. fasta sequence database)

3. Run the workflow:
    - To create a general spectral library: `snakemake . --use-conda`
    - To create a calibrated retention time dataset: `snakemake --snakefile make_rt_lib.smk --use-conda`

## Configuration
All settings can be found in a JSON configuration file: `conf/snakemake_config.json`.

| Section | Option | Default value | Description |
|---|---|---|---|
| download | pxd_identifier | "PXD000000" | PXD identifier of PRIDE Archive project to download. |
| | file_pattern | ".\*" | Regular expression that matches all raw file filenames to download (`.*` matches all filenames). |
| convert | exec | "ThermoRawFileParser.sh" | Executable command to call ThermoRawFileParser. See [Note 1](#note-1). |
| search | msgfplus_conf | "conf/msgfplus_params.txt" | Path to MSGFPlus configuration file. |
| | fasta | "path/to/search_db.fasta" | Path to protein fasta. Important: MSGFPlus will add decoy peptides by default; they should not yet be present in the given fasta file. |
| | msgfplus_exec | "msgf_plus" | Executable command to call MSGFPlus. See [Note 2](#note-2). |
| | threads_per_search | 5 | Number of threads per MSGFPlus search. See [Note 3](#note-3).

### Note 1
**ThermoRawFileParser executable**  
To use a custom version of the ThermoRawFileParser, change convert > exec to:
- Using environment TRFP: `"ThermoRawFileParser.sh"`
- Using custom TRFP: `"mono /path/to/ThermoRawFileParser.exe"`

### Note 2
**MSGFPlus executable**  
Idem for MSGFPlus, change search > msgfplus_exec to:
- Using environment MSGFPlus: `msgf_plus`
- Using custom jar file: `"msgfplus_exec": "java -Xmx5000M -jar /path/to/MSGFPlus/MSGFPlus.jar"`

The latter allows a custom memory limit for the Java VM. By default, this is 1GB.

### Note 3
**Indirectly limit memory usage while searching with threads_per_search**  
The option search > threads_per_search defines the number of threads each
individual search can use. In combination with the snakemake `--cores x` option,
this allows you to constrict the number of parallel searches. E.g.: the
combination of `--cores 24` and `threads_per_search: 6` limits the number of
parallel searches to 4. This can be convenient if you would run into memory
issues caused by too many parallel searches.
