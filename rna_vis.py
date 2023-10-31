from Bio.PDB import MMCIFParser
from structure_data import StructureData
from split_entities import splitEntities
import os, sys
from d3graph import d3graph, vec2adjmat
import pandas as pd
import json
import collections
from get_interactions import getInteractions
from run_dssr import runDSSR
from get_edges import getEdges
from process_graph import processEdges, processNodes
import json
from hbond_extractor import hbondExtractor, labelHbondEdges
import sys
from get_ss import getSS, processSS

parser = MMCIFParser()

#update: no need to change anymore
home =  os.path.dirname(os.path.abspath(__file__)) #change this line only 

pdb_path = "{}/dssr_output/".format(home)
# pdb_file = "8fvi-assembly1.cif"
prefix = '1ivs'

if len(sys.argv) > 1:
   prefix = sys.argv[1]

pdb_file = "{}.tmp.cif".format(prefix)
structure = StructureData(os.path.join(pdb_path, pdb_file), name="co_crystal")
protein, rna = splitEntities(structure) # split RNA and protein from structure

ss = getSS(prefix)
# print(ss)

with open("{}/{}-dssr.json".format(pdb_path, prefix)) as FH:
   data = json.load(FH, object_pairs_hook=collections.OrderedDict) 


protein_interactions,ss_dict = getInteractions(protein, rna, prefix)
pairs,backbone_edges, interaction_edges, interaction_types, stacks = getEdges(data, protein_interactions, ss_dict)

#update: added functions to extract all H-bond interactions from dssr and to add H-bond labels to interaction_types object
hbond_set = hbondExtractor(data)
interaction_types  = labelHbondEdges(interaction_types, hbond_set)

all_edges = pairs + backbone_edges + interaction_edges + stacks
# all_nodes = set()
# for tuple in all_edges:
#    if tuple[0] not in all_nodes:
#       all_nodes.add(tuple[0])
#    if tuple[1] not in all_nodes:
#       all_nodes.add(tuple[1])

# print(all_nodes)
# df_nodes = pd.DataFrame(all_nodes, columns=['id'])
# df_nodes_json = df_nodes.to_json(orient='records')
# print(df_nodes_json)

# df_links = pd.DataFrame(all_edges, columns=['source', 'target'])
# df_links['weight'] = [20]*len(pairs) + [100]*(len(backbone_edges)) + [5]*(len(interaction_edges)) + [20]*(len(stacks))
# df_links_json = df_links.to_json(orient='records')
# print(df_links_json)
d3 = d3graph(support=None, collision=0.5)
df = pd.DataFrame(all_edges, columns=['source', 'target'])
df['weight'] = [100]*len(pairs) + [100]*(len(backbone_edges)) + [5]*(len(interaction_edges)) + [20]*(len(stacks))
adjmat = vec2adjmat(df['source'], df['target'], weight=df['weight'])

###############################################################

d3.graph(adjmat)

d3.set_edge_properties(directed=True) # setting earlier to then update?


# can probably pre-compute, then add to the dataframe and use that?
# iterate through nodes to change colors, label, etc.

d3.node_properties = processNodes(d3.node_properties)
d3.edge_properties = processEdges(d3.edge_properties, backbone_edges, stacks, pairs, interaction_types)

# d3.show(filepath='{}/output/{}.html'.format(home, pdb_file), show_slider=False, showfig=False)
# click={'fill': None, 'stroke': '#F0F0F0', 'size': 2.5, 'stroke-width': 10} # add inside d3 show to highlight click
final_json = d3.show(filepath='{}/output/{}.html'.format(home, pdb_file), show_slider=False, showfig=False)
print(final_json)
final_json_str = json.dumps(final_json)
# print(final_json_str)

# Generate file for subgraph testing
with open('{}/output/{}_graph.json'.format(home, prefix), 'w') as outfile:
    json.dump(final_json, outfile)

ss_json = processSS(ss)
