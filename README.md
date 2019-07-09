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
- Conda
- Python 3
- Installed environment from `environment.yml`.

Create and activate the environment:
```
conda env create -f environment.yml
conda activate pxd_to_speclib
```

## Configuration
All settings can be found in a JSON configuration file: `conf/snakemake_config.json`.

To use a custom version of the ThermoRawFileParser, change convert > exec to:
- Using environment TRFP: `"ThermoRawFileParser.sh"`
- Using custom TRFP: `"mono /path/to/ThermoRawFileParser.exe"`

Idem for MSGFPlus, change search > msgfplus_exec to:
- Using environment MSGFPlus: `msgf_plus`
- Using custom jar file: `"msgfplus_exec": "java -Xmx5000M -jar /path/to/MSGFPlus/MSGFPlus.jar"`

The latter allows a custom memory limit for the Java VM. By default, this is 1GB.

The option search > threads_per_search defines the number of threads each
individual search can use. In combination with the snakemake `--cores x` option,
this allows you to constrict the number of parallel searches. E.g.: the
combination of `--cores 24` and `threads_per_search: 6` limits the number of
parallel searches to 4. This can be convenient if you would run into memory
issues caused by too many parallel searches.
