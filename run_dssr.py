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
    outpath = "/home/raktim/rnaprodb/rnaprodb/dssr_output/{}".format(tmpdir)
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    if not isinstance(structure, str):
        file_name = "{}/{}.tmp.cif".format(outpath, prefix)
        structure.save(file_name)
    else:
        file_name = structure
    
    args = ["x3dna-dssr", "--i={}".format(file_name), "--o={}/{}-dssr.json".format(outpath, prefix), "--json", "--more", "--idstr=long", "--non-pair"]
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
    
    fname = sys.argv[1]
    
    #update: no need to change anymore
    home =  os.path.dirname(os.path.abspath(__file__)) #change this line only 
    pdb_path = "/home/raktim/rnaprodb/data/cifs/"
    
    # pdb_file = "8fvi-assembly1.cif"
    
    for item in tqdm(open(fname,"r").readlines()):
        item = item.strip()
        #f = os.path.join(cifdir, fname)
        pdb_file = item #"{}-assembly1.cif".format(prefix)
        original_structure = StructureData(os.path.join(pdb_path, pdb_file), name="co_crystal")
        protein, rna = splitEntities(original_structure) # split RNA and protein from structure

        # NOTE: here we save the cleaned temporary file before processing it in runDSSR. For larger files, this poses a challenge.
        rna = cleanRNA(rna)


        full = rna[0]
        for chain in protein.get_chains():
            while chain.id in [i.id for i in full.child_list]:
                chain.id = chain.id+chain.id
            full.add(chain)
        full = StructureData(full)
        data = runDSSR(full, quiet=True, prefix=pdb_file.split("-")[0], tmpdir="")

