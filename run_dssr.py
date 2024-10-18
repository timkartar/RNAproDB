import subprocess
import json
import os, sys
import collections
from structure_data import StructureData
from split_entities import splitEntities
from clean_rna import cleanRNA
import os, sys
from tqdm import tqdm
from Bio.PDB import PDBIO, MMCIFIO, PDBParser, MMCIFParser

parser= MMCIFParser()

backend =  os.path.dirname(os.path.abspath(__file__))
cifs = backend + "/output/cifs/"

def runDSSR(structure, quiet=True, prefix='rna', tmpdir=''):
    """Run DSSR on given PDB file.
    
    Parameters
    ----------
    prefix: string
        The file prefix.
    """
    pdbpath = cifs + "/{}".format(tmpdir)
    outpath = backend + "/dssr_output/{}".format(tmpdir)
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    if not isinstance(structure, str):
        #print("prefix")
        #print(prefix)
        file_name = "{}/{}-assembly1.cif".format(pdbpath, prefix)
        #structure.save(file_name)
    else:
        file_name = structure
    
    args = ["x3dna-dssr", f"--i={file_name}", "--o={}/dssr_output/{}-dssr.json".format(backend,
        prefix), "--json", "--more", "--idstr=long", "--non-pair"]#, "--nested"]
    if quiet:
        FNULL = open(os.devnull, 'w')
        subprocess.call(args, stdout=FNULL, stderr=FNULL)
        subprocess.call(["x3dna-dssr", "--cleanup"],stdout=FNULL, stderr=FNULL)
        FNULL.close()
    else:
        subprocess.call(args)
        subprocess.call(["x3dna-dssr", "--cleanup"])
    
    with open("{}/{}-dssr.json".format(outpath, prefix)) as FH:
        DSSR = json.load(FH, object_pairs_hook=collections.OrderedDict)
    

    return DSSR

if __name__ == "__main__":

    pdb_id = sys.argv[1]
    
    #update: no need to change anymore
    # home =  os.path.dirname(os.path.abspath(__file__)) #change this line only
    pdb_path = cifs
    # pdb_file = f"{pdb_id}-assembly1.cif"
    
    for item in tqdm(open(cifs+"/{}-assembly1.cif".format(pdb_id),"r").readlines()):
        item = item.strip()
        #f = os.path.join(cifdir, fname)
        pdb_file = item #"{}-assembly1.cif".format(prefix)
        if os.path.exists(backend + "/dssr_output/{}-dssr.json".format(pdb_id)):
            print("here")
            print("{}-already exists".format(pdb_id))
            continue
        
        #try:
        print("there")
        original_structure = parser.get_structure(pdb_id,pdb_path)
        #original_structure = StructureData(os.path.join(pdb_path, f'{pdb_id}-assembly1.cif'), name="co_crystal")

        # NOTE: here we save the cleaned temporary file before processing it in runDSSR. For larger files, this poses a challenge.
        #rna = cleanRNA(rna)

        '''
        #protein, rna = splitEntities(original_structure) # split RNA and protein from structure
        full = rna[0]
        for chain in protein.get_chains():
            while chain.id in [i.id for i in full.child_list]:
                chain.id = chain.id+chain.id
            full.add(chain)
        full = StructureData(full)
        '''
        #data = runDSSR(full, quiet=True, prefix=pdb_id, tmpdir="")
        data = runDSSR(original_structure, quiet=True, prefix=pdb_id, tmpdir="")
        #except Exception as e:
        #    print(prefix, e)

