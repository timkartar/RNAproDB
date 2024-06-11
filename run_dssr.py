import subprocess
import json
import os, sys
import collections
from structure_data import StructureData
from split_entities import splitEntities
from clean_rna import cleanRNA
import os, sys
from tqdm import tqdm

def runDSSR(structure, quiet=True, prefix='rna', tmpdir=''):
    """Run DSSR on given PDB file.
    
    Parameters
    ----------
    prefix: string
        The file prefix.
    """
    outpath = "/mnt/c/Users/hosse/Desktop/Github/backend_master/dssr_output/{}".format(tmpdir)
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    if not isinstance(structure, str):
        print("prefix")
        print(prefix)
        file_name = "{}/{}.tmp.cif".format(outpath, prefix)
        structure.save(file_name)
    else:
        file_name = structure
    
    args = ["/mnt/c/Users/hosse/Desktop/Github/DeepPBS/dependencies/bin/x3dna-dssr", f"--i={file_name}", "--o=/mnt/c/Users/hosse/Desktop/Github/backend_master/dssr_output/{}-dssr.json".format(pdb_id), "--json", "--more", "--idstr=long", "--non-pair"]
    if quiet:
        FNULL = open(os.devnull, 'w')
        subprocess.call(args, stdout=FNULL, stderr=FNULL)
        subprocess.call(["/mnt/c/Users/hosse/Desktop/Github/DeepPBS/dependencies/bin/x3dna-dssr", "--cleanup"],stdout=FNULL, stderr=FNULL)
        FNULL.close()
    else:
        subprocess.call(args)
        subprocess.call(["/mnt/c/Users/hosse/Desktop/Github/DeepPBS/dependencies/bin/x3dna-dssr", "--cleanup"])
    
    with open("{}/{}-dssr.json".format(outpath, prefix)) as FH:
        DSSR = json.load(FH, object_pairs_hook=collections.OrderedDict)
    

    return DSSR

if __name__ == "__main__":

    pdb_id = sys.argv[1]
    
    #update: no need to change anymore
    # home =  os.path.dirname(os.path.abspath(__file__)) #change this line only 
    pdb_path = "/mnt/c/Users/hosse/Desktop/Github/frontend_master/public/cifs/"
    # pdb_file = f"{pdb_id}-assembly1.cif"
    
    for item in tqdm(open(f"/mnt/c/Users/hosse/Desktop/Github/frontend_master/public/cifs/{pdb_id}-assembly1.cif","r").readlines()):
        item = item.strip()
        #f = os.path.join(cifdir, fname)
        pdb_file = item #"{}-assembly1.cif".format(prefix)
        if os.path.exists("/mnt/c/Users/hosse/Desktop/Github/backend_master/dssr_output/{}-dssr.json".format(pdb_id)):
            print("here")
            print("{}-already exists".format(pdb_id))
            continue
        
        #try:
        print("there")
        original_structure = StructureData(os.path.join(pdb_path, f'{pdb_id}-assembly1.cif'), name="co_crystal")
        protein, rna = splitEntities(original_structure) # split RNA and protein from structure

        # NOTE: here we save the cleaned temporary file before processing it in runDSSR. For larger files, this poses a challenge.
        #rna = cleanRNA(rna)


        full = rna[0]
        for chain in protein.get_chains():
            while chain.id in [i.id for i in full.child_list]:
                chain.id = chain.id+chain.id
            full.add(chain)
        full = StructureData(full)
        data = runDSSR(full, quiet=True, prefix=pdb_id, tmpdir="")
        #except Exception as e:
        #    print(prefix, e)

