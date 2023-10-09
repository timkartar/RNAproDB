import subprocess
import json
import os, sys
import collections

def runDSSR(structure, quiet=True, prefix='rna'):
    """Run DSSR on given PDB file.
    
    Parameters
    ----------
    prefix: string
        The file prefix.
    """
    
    if not isinstance(structure, str):
        file_name = "./dssr_output/{}.tmp.mmcif".format(prefix)
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
    fname = sys.argv[1]
    cifdir = "../data/cifs/"
    for item in open(fname,"r").readlines():
        item = item.strip()
        f = os.path.join(cifdir, fname)

