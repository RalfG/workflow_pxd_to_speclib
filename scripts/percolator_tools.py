import re
import pandas as pd


def psmid_to_usi(psmid, project_id):
    """
    Convert Percolator out PSMId to HUPO-PSI Universal Spectrum Identifier.
    See http://www.psidev.info/usi for more info.
    
    Expects the following formatted PSMId:
    `run` _ `SII` _ `MSGFPlus spectrum index` _ `PSM rank` _ `scan number` _ `MSGFPlus-assigned charge` _ `rank`
    See https://github.com/percolator/percolator/issues/147
    """
    psmid = psmid.split('_')
    usi = ':'.join(['mzspec', project_id, '_'.join(psmid[:-6]), 'scan', psmid[-3]])    
    return usi


def psmid_to_run(psmid):
    """
    Extract run from Percolator PSMId.
    
    Expects the following formatted PSMId:
    `run` _ `SII` _ `MSGFPlus spectrum index` _ `PSM rank` _ `scan number` _ `MSGFPlus-assigned charge` _ `rank`
    See https://github.com/percolator/percolator/issues/147
    """ 
    psmid = psmid.split('_')
    run = '_'.join(psmid[:-6])
    return run


def psmid_to_charge(psmid):
    """
    Extract charge from Percolator PSMId.
    
    Expects the following formatted PSMId:
    `run` _ `SII` _ `MSGFPlus spectrum index` _ `PSM rank` _ `scan number` _ `MSGFPlus-assigned charge` _ `rank`
    See https://github.com/percolator/percolator/issues/147
    """ 
    psmid = psmid.split('_')
    charge = int(psmid[-2])
    return charge


def psmid_to_scan(psmid):
    """
    Extract charge from Percolator PSMId.
    
    Expects the following formatted PSMId:
    `run` _ `SII` _ `MSGFPlus spectrum index` _ `PSM rank` _ `scan number` _ `MSGFPlus-assigned charge` _ `rank`
    See https://github.com/percolator/percolator/issues/147
    """ 
    psmid = psmid.split('_')
    scan = int(psmid[-3])
    return scan


def fix_pin_tabs(path, prot_sep='|||'):
    """
    Take a pin file and rewrite it, replacing the tabs that separate the
    Proteins column with a different separator
    """
    f = open(path)
    rows = f.readlines()
    outfile = path + '_fixed'
    out = open(outfile, 'w+')

    for i, row in enumerate(rows):
        if i == 0 & row.startswith('SpecId'):
            numcol = len(row.split('\t'))
            out.write(row)
        elif i == 1 & row.startswith('DefaultDirection'):
            out.write(row)
        else:
            r = row.strip().split('\t')
            r_cols = r[:numcol-1]
            r_proteins = r[numcol-1:]
            r_cols.append(prot_sep.join(r_proteins))
            out.write('\t'.join(r_cols) + '\n')
    f.close()
    out.close()
    return None


def extract_seq_mods(df, mods):
    """
    Extract PEPREC-style modifications and sequence from Percolator-
    style peptide notation.
    """

    # Add modifications column to PEPREC file
    # the keys correspond to the UNIMOD keys for each modification
    modifications = {}
    for mod in mods:
        modifications[str(mod["unimod_accession"])] = mod["name"]

    modlist = []
    # TODO get rid of iterrows!
    for _, row in df.iterrows():
        if 'UNIMOD' in row['modified_peptide']:
            pep = row['modified_peptide'].split('.')[1]
            mods = re.findall(r'\[([^]]*)\]', pep)
            modstring = ''
            for mod in mods:
                mod = '[' + mod + ']'
                key = mod.split(':')[1].rstrip(']')
                try:
                    if key == '21':
                        phospholoc = pep[pep.find(mod)-1]
                        modstring += str(pep.find(mod)) + '|' + modifications[key] + phospholoc + '|'
                        pep = pep.replace(mod, '', 1)
                    else:
                        modstring += str(pep.find(mod)) + '|' + modifications[key] + '|'
                        pep = pep.replace(mod, '', 1)
                except:
                    print('Modification not expected: {}'.format(mod))
            modlist.append(modstring.rstrip('|'))
        else:
            modlist.append('')
    
    # Rewrite peptide sequences without the UNIMOD modification ids
    peplist = []
    for _, row in df.iterrows():
        pep = row['modified_peptide']
        pep = pep.split('.')[1]
        if 'UNIMOD' in pep:
            mods = re.findall(r'\[([^]]*)\]', pep)
            for mod in mods:
                pep = pep.replace('[' + mod + ']', '', 1)
        peplist.append(pep)

    df_out = pd.DataFrame({'peptide': peplist, 'modifications': modlist})
    return df_out
