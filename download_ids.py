#!/usr/bin/env python
import argparse

from os.path import join as ospj
import json
import requests
from string import Template


def main():
    ### Query RCSB for new entries
    search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
    query = open("./search_query.tmpl").read()
    query = Template(query)
    query = json.loads(query.substitute(release_date='1970-10-01'))

    req = requests.post(search_url, json=query)
    RESULTS = req.json()

    ### Download structures
    #root_download_path = ospj(ROOT_DIR, "external/PDB/pdb_entries")
    #download_url = CF["RCSB"]["mmcif_data_url"]
    #downloaded = open(os.path.join(ROOT_DIR, "external/PDB/pdb_entries", "newest_downloaded_releases.txt"), "w")
    for entry in RESULTS["result_set"]:
        entry_id  = entry["identifier"].lower()
        #path = os.path.join(root_download_path, entry_id[0], "{}.cif.gz".format(entry_id))
        print(entry_id)
        # download entry
    ## query no protein entries
    # query = open("./search_query_noprotein.tmpl").read()
    ## query all structures
    query = open("./search_query.tmpl").read()


    query = Template(query)
    query = json.loads(query.substitute(release_date='1970-10-01'))

    req = requests.post(search_url, json=query)
    RESULTS = req.json()

    ### Download structures
    #root_download_path = ospj(ROOT_DIR, "external/PDB/pdb_entries")
    #download_url = CF["RCSB"]["mmcif_data_url"]
    #downloaded = open(os.path.join(ROOT_DIR, "external/PDB/pdb_entries", "newest_downloaded_releases.txt"), "w")
    
    result_list =  []
    for entry in RESULTS["result_set"]:
        entry_id  = entry["identifier"].lower()
        #path = os.path.join(root_download_path, entry_id[0], "{}.cif.gz".format(entry_id))
        print(entry_id)
        # download entry
        result_list.append(entry_id)

    print(len(RESULTS["result_set"]))
    return result_list

if __name__ == "__main__":
    main()