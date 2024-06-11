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
    if(len(temp_split[1]) == 3):
        return True
    return False

"""
Returns ('p'/'n', name, position, chain, ss (protein only))
ASSUMES: proteins are represented by 3 letter code, and nucleotides by one letter!
"""
def parse_node(node_id):
    print(f"Parsing node_id: {node_id}")
    temp_split = node_id.split(':')
    print(temp_split)
    # Is a protein
    if (temp_split[1] in ['ZN','NA','MG']):
        return ('n', temp_split[1], temp_split[2], temp_split[0])
    elif(len(temp_split[1]) == 1 or (len(temp_split[1]) == 3 and temp_split[1] not in d3to1.keys())): #Is a nucleotide
        return ('n', temp_split[1], temp_split[2], temp_split[0])
    elif(is_a_protein(node_id)):
        #check if has secondary structure
        # if len(temp_split) == 3: # no secondary structure!
        #     # print(node_id)
        #     return ('p', temp_split[1], temp_split[2], temp_split[0], "Unknown")
        # else:
        return ('p', temp_split[1], temp_split[2], temp_split[0], temp_split[3])
    else: #ERROR
        print("Error parsing node!")
        return None
    
d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
        'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
        'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
        'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
