import os, sys
import subprocess
#from dnaprodb_utils import CHAIN_RE, RESN_RE, RESI_RE, ATOM_RE
#from dnaprodb_utils import residueMoiety
#from dnaprodb_utils import nucleotideMoiety
#from dnaprodb_utils import log, getHash, getID, getCM
import re
from Bio.PDB import PDBIO
import copy
from utilities import nt_colors, chem_components, d3to1
def is_nt(resname):
    if resname in nt_colors.keys() or resname in chem_components.keys():
        return True
    else:
        return False

def is_aa(resname):
    if resname in d3to1.keys():
        return True
    elif is_nt(resname):
        return False
    return 'x'

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
    rc = subprocess.call(['{}/external/hbplus'.format(backend), '-h', '3.0', '-d', '3.5', '{}.pdb'.format('/srv/www/rnaprodb/rnaprodb_dev/1ivs-assembly1'),
        '{}.pdb'.format('/srv/www/rnaprodb/rnaprodb_dev/1ivs-assembly1')])#, stdout=FNULL, stderr=FNULL)

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
        items["d_ins"] = d_ins.replace(" ","")
        items["d_atom"] = d_atom
        items["a_chain"] = a_chain
        items["a_chain"] = a_chain
        items["a_resi"] = a_resi
        items["a_resn"] = a_resn
        items["a_ins"] = a_ins.replace(" ","")
        items["a_atom"] = a_atom
        items["dist"] = dist

        if (d_resn == "HOH") != (a_resn == "HOH"):
            water_hbonds.append(items)
        else:
            pass

    os.remove('{}.pdb'.format(prefix))
    os.remove('{}.hb2'.format(prefix))
    return water_hbonds

def getWHbonds(path, prefix, structure, ss_dict, interaction_types):
    data = runHBplus(path, prefix, structure)
    waters = {}
    for item in data:
        other_data = {}
        water = "{}:{}:{}".format(item["a_chain"], item["a_resi"], item["a_ins"])
        other = "{}:{}:{}".format(item["d_chain"], item["d_resi"], item["d_ins"])
        if item["d_resn"] == "HOH":
            water, other = other, water
            othername = item["a_resn"] 
            other_data['type'] = 'acceptor'
            other_data['atom'] = item['a_atom']
        else:
            othername = item["d_resn"] 
            other_data['type'] = 'donor'
            other_data['atom'] = item['d_atom']

        other_data['distance_to_water'] = item['dist']
        other_data['rnaprodb_id'] = other
        other_data['name'] = othername
    
        if water not in waters.keys():
            waters[water] = [[],[]]
        if is_aa(othername) == 'x': ##ignore modified aa e.g. MSE
            #print(othername, is_aa(othername))
            continue
        if is_aa(othername):
            other_data['rnaprodb_id'] = ":".join(other.split(":")[:-1]) + ":"
            #waters[water][1].append(other)
            cand = [other, other_data]
            #if cand not in waters[water][1]:
            waters[water][1].append(cand)
        elif is_nt(othername):
            cand = [other, other_data]
            #if cand not in waters[water][0]:
            waters[water][0].append(cand)
            #waters[water][0].append(other)
    #exit()
    keys = list(waters.keys())
    for water in keys:
        if ([] in waters[water]):
            del waters[water]

    wm_edges = set()
    for item in waters.values():
        for nt_data in item[0]:
            nt_data = nt_data[1]
            ntsplit = nt_data['rnaprodb_id'].split(":")
            for aa_data in item[1]:
                aa_data = aa_data[1]
                aasplit = aa_data['rnaprodb_id'].split(":")
                nt_node = ":".join([ntsplit[0], nt_data['name'], ntsplit[1], ntsplit[2]])
                aa_node = ":".join([aasplit[0], aa_data['name'], aasplit[1],
                        ss_dict[aa_data['rnaprodb_id']]])
                edge = (nt_node, aa_node)
                wm_edges.add(edge)
                if edge not in interaction_types:
                    interaction_types[edge] = {'whbond'}
                else:
                    interaction_types[edge].add('whbond')
    wm_edges = list(wm_edges)
    print(wm_edges)
    return wm_edges, interaction_types, waters
if __name__=="__main__":
    pass
    #runHBplus()
