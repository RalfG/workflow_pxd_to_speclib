"""
## Download PRIDE Project
Download PRIDE project files for a given PRIDE identifier. With the `-f`
argument certain file types can be chosen for download.

*Download_PRIDE_Project.py*
**Input:** PRIDE Archive identifier
**Output:** Downloaded files, sorted in folders by file type
"""

import os
import json
import requests
import argparse

import wget
import pandas as pd


def argument_parser():
    parser = argparse.ArgumentParser(description='Download files from PRIDE Archive for a given project.')
    parser.add_argument('pxd_identifier', action='store',
                        help='PXD identifier of project from which to download files')
    parser.add_argument('-p', dest='pattern', action='store',
                        help='Rexex pattern matching to files to be downloaded')
    parser.add_argument('-f', dest='filetypes', action='store', nargs='+',
                        help='filetypes to download (msf, raw, txt, zip...)')
    args = parser.parse_args()
    return args


def get_files_df(pxd_identifier, filetypes, pattern):
    """
    Get DataFrame with files to download, filtered by extension and filename
    regex pattern.
    """
    get_files_url = "https://www.ebi.ac.uk:443/pride/ws/archive/file/list/project"

    # Get dataframe with file info through PRIDE Archive REST API
    url = "{}/{}".format(get_files_url, pxd_identifier)
    response = pd.DataFrame(json.loads(requests.get(url).content.decode('utf-8'))['list'])
    response['fileExtension'] = response['fileName'].str.split('.').apply(lambda x: x[-1])

    # Set regex pattern
    if pattern:
        response = response[response['fileName'].str.contains(pattern)]

    # Set extensions
    extensions = response['fileExtension'].unique()
    #print("Found files with extensions: \t\t{}".format(extensions))
    if filetypes:
        extensions = filetypes
        response = response[response['fileExtension'].str.lower().isin([x.lower() for x in filetypes])]

    #print("Downloading files with extensions: \t{}".format(extensions))

    return response


def run():
    args = argument_parser()
    get_meta_url = "https://www.ebi.ac.uk:443/pride/ws/archive/project"

    # Make folder for project and download meta data
    print("Downloading meta data...")
    url = "{}/{}".format(get_meta_url, args.pxd_identifier)
    response = json.loads(requests.get(url).content.decode('utf-8'))
    with open("pxd_project_metadata.json", "w") as f:
        f.write(json.dumps(response))
    with open("pxd_project_metadata.txt", "w") as f:
        f.write('{}\n'.format(response['accession']))
        f.write('{}\n'.format(response['publicationDate']))
        f.write('{}\n'.format(response['title']))
        try:
            f.write('{}\n'.format(response['references'][0]['desc']))
            for _, item in enumerate(response['references'][0]['ids']):
                f.write(item + '\n')
        except IndexError:
            pass

    # Download files
    print("Downloading files...")
    response = get_files_df(args.pxd_identifier, args.filetypes, args.pattern)
    for ext in args.filetypes:
        if not os.path.exists(ext):
            os.mkdir(ext)
        count = 0
        total = len(response[response['fileExtension'] == ext])
        for _, row in response[response['fileExtension'] == ext].iterrows():
            count += 1
            #if count % 10 == 0 or count == 1:
            print(" | {} {}/{}".format(ext, count, total), end=' ')
            target_path = os.path.join(ext, row['fileName'])
            if not os.path.isfile(target_path):
                wget.download(row['downloadLink'], out=ext)
                #os.system("wget -O '{}' '{}' -q --show-progress".format(target_path, row['downloadLink']))


if __name__ == '__main__':
    run()
