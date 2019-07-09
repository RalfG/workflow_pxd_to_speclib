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

To use a custom version of the ThermoRawFileParser, change
```
"convert": {
    "exec": "ThermoRawFileParser.sh"
}
```
to
```
"convert": {
    "exec": "mono /path/to/ThermoRawFileParser.exe"
}`
```
in the JSON config file.
