from Bio.PDB import is_aa
import numpy as np
import os
chem_components = dict(np.load(os.path.dirname(os.path.abspath(__file__)) + "/modified_parents.npz",allow_pickle=True))
d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
        'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
        'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
        'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
nt_colors = {'A': '#FF9896',#'#90cc84',
    'C': '#DBDB8D',#'#AEC7E8',
    'G': '#90cc84',#'#DBDB8D',
    'U': '#AEC7E8',#'#FF9896',
    'T': '#AEC7E8',#'#FF9896',
    'DA': '#FF9896',#'#90cc84',
    'DC': '#DBDB8D',#'#AEC7E8',
    'DG': '#90cc84',#'#DBDB8D',
    'DT': '#AEC7E8'#'#FF9896'

}

"""
Input: ('p'/'nt', name, position, chain, ss (protein only))
Returns text string from node tuple representation
"""
def node_to_text(node):
    if len(node) == 5: # IS A PROTEIN
        return node[3] + ":" + node[1] + ":" + node[2] + ":" + node[4]
    return node[3] + ":" + node[1] + ":" + node[2] # is a NT

"""
Returns rnaprodb nucleotide text string (A:B:C) from dssr id (1..C.C.955.) representation
"""
def dssr_id_to_text(dssr_id):
    return ":".join(dssr_id.split(".")[2:-1])

"""
Checks if node_id (A:B:C) is a protein
"""
def is_a_protein(node_id):
    temp_split = node_id.split(':')
    if(temp_split[1] in d3to1.keys()):
        return True
    return False

"""
Returns ('p'/'n', name, position, chain, ss (protein only))
ASSUMES: proteins are represented by 3 letter code, and nucleotides by one letter!
"""
def parse_node(node_id):
    #print(f"Parsing node_id: {node_id}")
    temp_split = node_id.split(':')
    #print(temp_split)
    # Is a protein
    if (temp_split[1] in ['ZN','NA','MG']):
        return ('', temp_split[1], temp_split[2], temp_split[0])
    elif(temp_split[1] in nt_colors.keys() or temp_split[1] in chem_components.keys()): #Is a nucleotide
        if temp_split[1] == "5MU":
            print(('n', temp_split[1], temp_split[2], temp_split[0]))
        return ('n', temp_split[1], temp_split[2], temp_split[0])
    elif(is_a_protein(node_id)):
        #check if has secondary structure
        # if len(temp_split) == 3: # no secondary structure!
        #     # print(node_id)
        #     return ('p', temp_split[1], temp_split[2], temp_split[0], "Unknown")
        # else:
        return ('p', temp_split[1], temp_split[2], temp_split[0], temp_split[3])
    else: #ERROR
        return ('x', temp_split[1], temp_split[2], temp_split[0])
        print("Error parsing node!")
        return None
    

