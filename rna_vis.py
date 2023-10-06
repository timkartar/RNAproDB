from Bio.PDB import MMCIFParser
from structure_data import StructureData
from split_entities import splitEntities
from clean_rna import cleanRNA
import os, sys
from run_dssr import runDSSR
from d3blocks import D3Blocks
import pandas as pd
import json
import collections
from get_interactions import getInteractions

nt_colors = {'A': '#00994C',
    'C': '#000099',
    'G': '#FBA922',
    'U': '#990000'
}

"""
Returns rnaprodb text string (A:B:C) from dssr id representation
"""
def dssr_id_to_text(dssr_id):
    return ":".join(dssr_id.split(".")[2:-1])

"""
Returns text string from node tuple representation
"""
def node_to_text(node):
    return node[3] + ":" + node[1] + ":" + node[2]

"""
IP: Returns if two nodes (tuple representation) are a backbone
"""
def is_backbone_edge(node_1, node_2):
    #if (both nucleotides) AND (one nt away from each other) AND (on same strand)
    if(node_1[0] == 'nt' and node_2[0] == 'nt' and ((int(node_2[2]) - int(node_1[2]) == 1) or (int(node_2[2]) - int(node_1[2]) == -1)) and node_1[3] == node_2[3]):
       return True
    else:
        return False
    
"""
stack[nts_long]: '..C.G.901.,..C.A.972.'
"""
def get_stacking_interactions(dssr):
    stack_interactions = []
    for stack in dssr['stacks']:
        nts_long_split = stack['nts_long'].split(',')
        first_nucleotide = dssr_id_to_text(nts_long_split[0])
        sec_nucleotide = dssr_id_to_text(nts_long_split[1])
        stack_interactions.append((first_nucleotide, sec_nucleotide))
    return stack_interactions


"""
Returns ('p'/'nt', name, position, chain)
ASSUMES: proteins are represented by 3 letter code, and nucleotides by one letter!
"""
def parse_node(node_id):
    temp_split = node_id.split(':')
    # Is a protein
    if(len(temp_split[1]) == 3):
        return ('p', temp_split[1], temp_split[2], temp_split[0])
    elif(len(temp_split[1]) == 1): #Is a nucleotide
        return ('nt', temp_split[1], temp_split[2], temp_split[0])
    else: #ERROR
        print("Error parsing node!")
        return None    

"""
Returns (('p'/'nt', name, position, chain),('p'/'nt', name, position, chain))
ASSUMES: proteins are represented by 3 letter code, and nucleotides by one letter!
"""
def parse_edge(node_id):
    first_node = parse_node(node_id[0])
    sec_node = parse_node(node_id[1])
    return first_node,sec_node

"""
Returns base and backbone pairings from pre-processed (i.e., RNA only) DSSR 
"""
def get_edges(dssr, protein_interactions):
    pairs_dict = {}
    pairs=[]

    #add pairs to dictionary for easy lookup
    for pair in data["pairs"]:
        p1 = dssr_id_to_text(pair.get("nt1"))
        p2 = dssr_id_to_text(pair.get("nt2"))
        pairs_dict[p1] = p2
        pairs_dict[p2] = p1
    # add base pairing and self edges
    
    for nt in data["nts"]:
        nt_id = dssr_id_to_text(nt.get("nt_id"))
        if nt_id in pairs_dict:
            # base pairing edge
            pair_edge = (nt_id, pairs_dict[nt_id])
        else:
            # self edge
            pair_edge = (nt_id, nt_id)
        pairs.append(pair_edge)
    # add backbone edges
    backbone_edges = []
    for i, nt in enumerate(data["nts"]):
        nt_id = dssr_id_to_text(nt.get("nt_id"))
        try:
            nt_next = dssr_id_to_text(data["nts"][i+1].get("nt_id"))
        except Exception as e:
            print(nt_id, e)
            continue
        c1 = nt_id.split(":")[0]
        c2 = nt_next.split(":")[0]

        if c1 == c2:
            backbone_edges.append((nt_id, nt_next))
        
    interaction_edges = []
    for key, val in protein_interactions.items():
        for v in val:
            nt = ":".join((key.split(":")[:-1]))
            interaction_edges.append((nt, v))
            
    interaction_edges = list(set(interaction_edges))
    stacks = get_stacking_interactions(dssr)

    return pairs,backbone_edges, interaction_edges,stacks

parser = MMCIFParser()

home = "/home/aricohen/Desktop/rnaprodb_dev/" #change this line only

pdb_path = "{}/pdb/".format(home)
# pdb_file = "8fvi-assembly1.cif"
pdb_file = "1ivs-assembly1.cif"
original_structure = StructureData(os.path.join(pdb_path, pdb_file), name="co_crystal")
protein, rna = splitEntities(original_structure) # split RNA and protein from structure

rna = cleanRNA(rna)
protein = protein #no need for cleanProtein at the moment
data = runDSSR(rna, quiet=True, prefix='1ivs')


## protein_dssr = runDSSR(prefix="1ivs_combined")
## get secondary structure for protein residues
## for protein interaction edges, color the protein nodes by their secondary structure ()

protein_interactions = getInteractions(protein, rna)


pairs,backbone_edges, interaction_edges,stacks = get_edges(data, protein_interactions)


all_edges = pairs + backbone_edges + interaction_edges + stacks

df = pd.DataFrame(all_edges, columns=['source', 'target'])
df['weight'] = [20]*len(pairs) + [100]*(len(backbone_edges)) + [5]*(len(interaction_edges)) + [20]*(len(stacks))
d3 = D3Blocks()
d3.d3graph(df, filepath='./')
d3.D3graph.set_edge_properties(directed=True, marker_color='red') # setting earlier to then update?

# can probably pre-compute, then add to the dataframe and use that?
# iterate through nodes to change colors, label, etc.
for node in d3.D3graph.node_properties:
    parsed_node = parse_node(node) # ('p'/'nt', name, position, chain)
    
    name = parsed_node[1]
    pos = str(parsed_node[2])
    chain = parsed_node[3]
    
    # global changes
    d3.D3graph.node_properties[node]['size'] = 12 # original 20
    d3.D3graph.node_properties[node]['label']= "---" + name # empty label, use tooltip instead
    d3.D3graph.node_properties[node]['opacity']= 0.705
    d3.D3graph.node_properties[node]['fontsize']= 20
    d3.D3graph.node_properties[node]['edge_size']= 1 # original 5

    if(parsed_node[0] == 'nt'): # is a nucleotide
        d3.D3graph.node_properties[node]['color']= nt_colors[name] #use nt color scheme
        tooltip = 'Nucleotide: ' + name +"\nPosition: " + pos + "\nChain: " + parsed_node[3]
        d3.D3graph.node_properties[node]['fontcolor']= nt_colors[name]
    else: # is protein residue
        d3.D3graph.node_properties[node]['size'] = 5 # original 5
        d3.D3graph.node_properties[node]['color']= 'gray' #use gray
        d3.D3graph.node_properties[node]['label']= '' # empty label, use tooltip instead
        tooltip = 'Resiude: ' + name +"\nPosition: " + pos + "\nChain: " + parsed_node[3]
    d3.D3graph.node_properties[node]['tooltip']= tooltip

# iterate through edges to determine colors, backbone edges 
for edge in d3.D3graph.edge_properties:
    first_node,sec_node = parse_edge(edge)
    edge_tuple = (node_to_text(first_node),node_to_text(sec_node)) #turn back into text to compare to backbone edge
    if edge_tuple in backbone_edges: #is a backbone edge. NOTE change to iterate through backbone edges instead!
        d3.D3graph.edge_properties[edge]['marker_start'] = ''
        # d3.D3graph.edge_properties[edge]['marker_end'] = 'arrow' # already set in set edge properties
        d3.D3graph.edge_properties[edge]['color'] = 'red' # works!
        # d3.D3graph.edge_properties[edge]['label_color'] = 'red'
        # d3.D3graph.edge_properties[edge]['label_fontsize'] = 8
        # d3.D3graph.edge_properties[edge]['marker_color'] = 'red' # BROKEN!
    else:
        d3.D3graph.edge_properties[edge]['marker_end'] = ''
        # d3.D3graph.edge_properties[edge]['directed'] = False

print(d3.D3graph.edge_properties)
d3.D3graph.show(filepath='{}/output/{}.html'.format(home, pdb_file))

# 'weight': 100.0, 'weight_scaled': 2.8129, 'edge_distance': 91.3043, 'color': '#808080', 'marker_start': '', 'marker_end': 'arrow', 'marker_color': '#808080', 'label': '', 'label_color': '#808080', 'label_fontsize': 8}