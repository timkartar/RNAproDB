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

nt_colors = {'A': '#00994C',
    'C': '#000099',
    'G': '#FBA922',
    'U': '#990000'
}

"""
Returns (nucleotide name, position)
"""
def parse_ntid(ntid):
    temp_split = ntid.split('.')
    return temp_split[3], str(temp_split[4])
    

"""
Returns base and backbone pairings from pre-processed (i.e., RNA only) DSSR 
"""
def get_edges(dssr, protein_interactions):
    pairs_dict = {}
    pairs=[]

    #add pairs to dictionary for easy lookup
    for pair in data["pairs"]:
        p1 = ":".join(pair.get("nt1").split(".")[2:-1])
        p2 = ":".join(pair.get("nt2").split(".")[2:-1])
        pairs_dict[p1] = p2
        pairs_dict[p2] = p1
    # add base pairing and self edges
    for nt in data["nts"]:
        nt_id = ":".join(nt.get("nt_id").split(".")[2:-1])
        
        if nt_id in pairs_dict:
            # base pairing edge
            pair_edge = (nt_id, pairs_dict[nt_id])
        else:
            # self edge
            pair_edge = (nt_id, nt_id)
        pairs.append(pair_edge)

    # add backbone edges
    backbone_edges = []
    for i in range(len(pairs)-1):
        backbone_edges.append((pairs[i][0], pairs[i+1][0]))
    
    
    interaction_edges = []
    for key, val in protein_interactions.items():
        for v in val:
            nt = ":".join((key.split(":")[:-1]))
            interaction_edges.append((nt, v))
            

    return pairs,backbone_edges, interaction_edges

# parser = MMCIFParser()

# home = "/Users/aricohen/Desktop/"
# pdb_path = "{}/rnaprodb_dev/".format(home)
# pdb_file = "1ivs-assembly1.cif"
# original_structure = StructureData(os.path.join(pdb_path, pdb_file), name="co_crystal")
# protein, rna = splitEntities(original_structure) # split RNA and protein from structure

rna = cleanRNA(rna)
protein = protein #no need for cleanProtein at the moment
data = runDSSR(rna, quiet=True, prefix='1ivs')
protein_interactions = getInteractions(protein, rna)

pairs,backbone_edges, interaction_edges = get_edges(data, protein_interactions)


all_edges = pairs + backbone_edges + interaction_edges

df = pd.DataFrame(all_edges, columns=['source', 'target'])
df['weight'] = [20]*len(pairs) + [40]*(len(backbone_edges)) + [5]*(len(interaction_edges))
d3 = D3Blocks()
d3.d3graph(df, filepath='./')

# iterate through nodes to change colors, label, etc.
for node in d3.D3graph.node_properties:
    nt_tuple = parse_ntid(node) # (nt, position)
    nt_name = nt_tuple[0]
    nt_pos = str(nt_tuple[1])

    d3.D3graph.node_properties[node]['size'] = 20
    d3.D3graph.node_properties[node]['label']= nt_pos
    d3.D3graph.node_properties[node]['color']= nt_colors[nt_name]
    d3.D3graph.node_properties[node]['opacity']= 0.705
    d3.D3graph.node_properties[node]['fontcolor']= '#000000'
    d3.D3graph.node_properties[node]['fontsize']= 20
    d3.D3graph.node_properties[node]['edge_size']= 5

d3.D3graph.show()