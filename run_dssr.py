import subprocess
import json
import os, sys
import collections
from Bio.PDB import MMCIFParser
from structure_data import StructureData
from split_entities import splitEntities
from clean_rna import cleanRNA
import os, sys

def runDSSR(structure, quiet=True, prefix='rna'):
    """Run DSSR on given PDB file.
    
    Parameters
    ----------
    prefix: string
        The file prefix.
    """
    
    if not isinstance(structure, str):
        file_name = "./dssr_output/{}.tmp.cif".format(prefix)
        structure.save(file_name)
    else:
        file_name = structure
    
    args = ["x3dna-dssr", "--i={}".format(file_name), "--o=./dssr_output/{}-dssr.json".format(prefix), "--json", "--more", "--idstr=long", "--non-pair"]
    if quiet:
        FNULL = open(os.devnull, 'w')
        subprocess.call(args, stdout=FNULL, stderr=FNULL)
        subprocess.call(["x3dna-dssr", "--cleanup"],stdout=FNULL, stderr=FNULL)
        FNULL.close()
    else:
        subprocess.call(args)
        subprocess.call(["x3dna-dssr", "--cleanup"])
    
    with open("./dssr_output/{}-dssr.json".format(prefix)) as FH:
        DSSR = json.load(FH, object_pairs_hook=collections.OrderedDict)
    

    return DSSR

if __name__ == "__main__":
    '''
    fname = sys.argv[1]
    cifdir = "../data/cifs/"
    parser = MMCIFParser()
    '''
    #update: no need to change anymore
    home =  os.path.dirname(os.path.abspath(__file__)) #change this line only 
    
    pdb_path = "{}/pdb/".format(home)
    # pdb_file = "8fvi-assembly1.cif"
    prefix = '1ivs'
    pdb_file = "{}-assembly1.cif".format(prefix)
    original_structure = StructureData(os.path.join(pdb_path, pdb_file), name="co_crystal")
    protein, rna = splitEntities(original_structure) # split RNA and protein from structure

    # NOTE: here we save the cleaned temporary file before processing it in runDSSR. For larger files, this poses a challenge.
    rna = cleanRNA(rna)

    '''
    for item in open(fname,"r").readlines():
        item = item.strip()
        f = os.path.join(cifdir, fname)
    '''
    full = rna[0]
    for chain in protein.get_chains():
        full.add(chain)
    
    data = runDSSR(rna, quiet=True, prefix='1ivs')

