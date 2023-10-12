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

parser = MMCIFParser()

#update: no need to change anymore
home =  os.path.dirname(os.path.abspath(__file__)) #change this line only 

pdb_path = "{}/dssr_output/".format(home)
# pdb_file = "8fvi-assembly1.cif"
prefix = '1ivs'
pdb_file = "{}.tmp.cif".format(prefix)
structure = StructureData(os.path.join(pdb_path, pdb_file), name="co_crystal")
protein, rna = splitEntities(structure) # split RNA and protein from structure

#   I updated the code to move running dssr to run_dssr and I run it now together on 
#   protein and RNA, this probably results in some changes in the output which breaks code
#   downstream in this file. 
#   TODO: uncomment the two lines below this comment, comment the runDSSR line and fix the 
#   breakage that happens then. Remove runDSSR import when it's done.

with open("{}/{}-dssr.json".format(pdb_path, prefix)) as FH:
   data = json.load(FH, object_pairs_hook=collections.OrderedDict) #TODO: uncomment this

# data = runDSSR(rna, quiet=True, prefix='1ivs') #TODO: comment this out

protein_interactions,ss_dict = getInteractions(protein, rna, prefix)
pairs,backbone_edges, interaction_edges,stacks = getEdges(data, protein_interactions, ss_dict)
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
d3 = d3graph()
df = pd.DataFrame(all_edges, columns=['source', 'target'])
df['weight'] = [20]*len(pairs) + [100]*(len(backbone_edges)) + [5]*(len(interaction_edges)) + [20]*(len(stacks))
adjmat = vec2adjmat(df['source'], df['target'], weight=df['weight'])
d3.graph(adjmat)
d3.set_edge_properties(directed=True, marker_color='red') # setting earlier to then update?


# can probably pre-compute, then add to the dataframe and use that?
# iterate through nodes to change colors, label, etc.

d3.node_properties = processNodes(d3.node_properties)
d3.edge_properties = processEdges(d3.edge_properties, backbone_edges, stacks)

d3.show(filepath='{}/output/{}.html'.format(home, pdb_file), click={'fill': None, 'stroke': '#F0F0F0', 'size': 2.5, 'stroke-width': 10})