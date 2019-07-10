"""
Percolator out (pout) to spectral library

Parse Percolator out and corresponding MGF files and filter to a spectral
library with unique peptides (seq, mods, charge).

JSON mods_config_file example:
```
{
    "modifications": [
        {"name":"Acetyl", "unimod_accession":1},
        {"name":"Oxidation", "unimod_accession":35},
        {"name":"Carbamidomethyl", "unimod_accession":4}
    ]
}
```
"""

# Standard library
import os
import json
import argparse
from glob import glob

# Third party
import pandas as pd

# Project
import percolator_tools
from parse_mgf import parse_mgf


def argument_parser():
    parser = argparse.ArgumentParser(description='Generate spectral library\
        from Percolator out (pout) and MGF files.')
    parser.add_argument('-c', dest='mods_config_file', action='store',
                        help='Path to JSON file with modifications info.',
                        required=True)
    parser.add_argument('-i', dest='project_id', action='store',
                        help='Project identifier (e.g. PXD). Required for\
                        Universal Spectrum Identifier.',
                        required=True)
    parser.add_argument('-p', dest='pout_path', action='store',
                        help='Path to directory with pout files.',
                        required=True)
    parser.add_argument('-m', dest='mgf_path', action='store',
                        help='Path to directory with MGF files.',
                        required=True)
    parser.add_argument('-o', dest='output_path', action='store',
                        help='Path to directory to write output files.',
                        required=True)
    parser.add_argument('-t', dest='fdr_threshold', action='store',
                        default=0.01, type=float, 
                        help='FDR treshold. PSMs with a q-value higher than the\
                        threshold are excluded from the spectral library.')
    parser.add_argument('-a', dest='all_spectra', action='store_true',
                        help='Do not filter for unique peptides (sequence,\
                        charge, modifications): include all spectra.')
    args = parser.parse_args()

    return args


def main():
    args = argument_parser()

    # Read JSON file with modifications:
    with open(args.mods_config_file) as json_file:  
        mods = json.load(json_file)['modifications']
    
    # Fix protein column in pout files
    for f in glob(os.path.join(args.pout_path, '*.pout')):
        percolator_tools.fix_pin_tabs(f)

    # Load all pout_fixed files
    all_pout_f = glob(os.path.join(args.pout_path, '*.pout_fixed'))
    to_concat = []
    for f in all_pout_f:
        df = pd.read_csv(f, sep='\t')
        #df['run'] = os.path.basename(f).split('.')[0]
        to_concat.append(df)
    all_pout = pd.concat(to_concat, axis=0)

    col_rename = {'peptide': 'modified_peptide', 'proteinIds': 'proteins', 'PSMId': 'percolator_psmid'}
    all_pout = all_pout.rename(columns=col_rename)

    # Parse required columns
    all_pout['usi'] = all_pout['percolator_psmid'].apply(lambda x: percolator_tools.psmid_to_usi(x, args.project_id))
    all_pout['run'] = all_pout['percolator_psmid'].apply(percolator_tools.psmid_to_run)
    all_pout['scan_number'] = all_pout['percolator_psmid'].apply(percolator_tools.psmid_to_scan)
    all_pout['charge'] = all_pout['percolator_psmid'].apply(percolator_tools.psmid_to_charge)

    # Filter on FDR threshold
    all_pout = all_pout[all_pout['q-value'] < args.fdr_threshold].copy()

    # Filter for best spectrum per peptide
    if not args.all_spectra:
        all_pout = all_pout.sort_values('q-value', ascending=True)
        all_pout = all_pout[~all_pout.duplicated(['modified_peptide', 'charge'], keep='first')].copy()
        all_pout = all_pout.sort_index().reset_index(drop=True)

    # Extract peptide and modifications out of modified_peptide column
    peprec_cols = percolator_tools.extract_seq_mods(all_pout, mods)
    all_pout = pd.concat([all_pout, peprec_cols], axis=1)

    # Parse all MGF files into one MGF with selected spectra
    parse_mgf(all_pout, args.mgf_path, outname=os.path.join(args.output_path + '/spectral_library.mgf'),
              filename_col='run', spec_title_col='scan_number',
              title_parsing_method='scan=', new_title_col='usi',
              show_progress_bar=False)

    # Create MS2PIP PEPREC (peptide record)
    peprec_cols = [
        'usi', 'modifications', 'peptide', 'charge',
        'proteins', 'score', 'q-value',
        'posterior_error_prob', 'run', 'scan_number'
    ] 
    peprec = all_pout[peprec_cols].rename(columns={'usi': 'spec_id'}).sort_values('scan_number')
    peprec.to_csv(os.path.join(args.output_path + '/spectral_library.peprec'), sep=' ', index=False)


if __name__ == '__main__':
    main()
