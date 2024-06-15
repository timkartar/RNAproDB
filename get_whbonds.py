import os, sys
import subprocess
#from dnaprodb_utils import CHAIN_RE, RESN_RE, RESI_RE, ATOM_RE
#from dnaprodb_utils import residueMoiety
#from dnaprodb_utils import nucleotideMoiety
#from dnaprodb_utils import log, getHash, getID, getCM
import re
from Bio.PDB import PDBIO
import copy
io = PDBIO()
backend =  os.path.dirname(os.path.abspath(__file__))
def runHBplus(path, prefix, structure):
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
    chain_map = {}
    structure = copy.deepcopy(structure)
    for chain in structure.get_chains():
        if len(chain.id) == 1:
            chain_map[chain.id] = chain.id
            alphabet = alphabet.replace(chain.id,"")
        else:
            chain_map[alphabet[0]] = chain.id
            chain.id = alphabet[0]
            alphabet = alphabet.replace(chain.id,"")



    io.set_structure(structure)
    io.save('{}.pdb'.format(prefix))
    rc = subprocess.call(['{}/external/hbplus'.format(backend), '-h', '3.0', '-d', '3.5', '{}.pdb'.format(prefix),
        '{}.pdb'.format(prefix)])#, stdout=FNULL, stderr=FNULL)

    water_hbonds = []
    HB = open('{}.hb2'.format(prefix),'r').readlines()
    for i in range(8,len(HB)):
        d_chain = chain_map[HB[i][0]]
        d_resi = str(int(HB[i][1:5].strip()))
        d_resn = HB[i][6:9].strip()
        d_ins = HB[i][5].replace('-',' ')
        d_atom = HB[i][9:13].strip()
        a_chain = chain_map[HB[i][14]]
        a_resi = str(int(HB[i][15:19].strip()))
        a_ins = HB[i][19].replace('-',' ')
        a_resn = HB[i][20:23].strip()
        a_atom = HB[i][23:27].strip()
        dist = float(HB[i][27:32].strip())

        items = {}
        items["d_chain"] = d_chain
        items["d_resi"] = d_resi
        items["d_resn"] = d_resn
        items["d_ins"] = d_ins
        items["d_atom"] = d_atom
        items["a_chain"] = a_chain
        items["a_chain"] = a_chain
        items["a_resi"] = a_resi
        items["a_resn"] = a_resn
        items["a_ins"] = a_ins
        items["a_atom"] = a_atom
        items["dist"] = dist

        if (d_resn == "HOH") != (a_resn == "HOH"):
            water_hbonds.append(items)
        else:
            pass

    os.remove('{}.pdb'.format(prefix))
    os.remove('{}.hb2'.format(prefix))
    return water_hbonds
if __name__=="__main__":
    pass
    #runHBplus()
