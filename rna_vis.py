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
from get_pca import getChainsAndPca, addPcaToGraph
from get_rnascape import addRNAscapeToGraph
from get_viennarna import addViennaToGraph
from get_num_nucleotides import count_nucleotides_slow, count_nucleotides_fast
from get_lw import getLW

parser = MMCIFParser()
home =  os.path.dirname(os.path.abspath(__file__))

pdb_path = "{}/dssr_output/".format(home)
# pdb_file = "8fvi-assembly1.cif"
prefix = '7vnv'

if len(sys.argv) > 1:
   prefix = sys.argv[1]

pdb_file = "{}.tmp.cif".format(prefix)
TOO_LARGE = False

# is too large flag. If too large, this script still pre-processes, but will flag it for views.py to know that a subgraph must be selected.
num_nts = count_nucleotides_fast(os.path.join(pdb_path, pdb_file))
if(num_nts > 500):
   TOO_LARGE = True

structure = StructureData(os.path.join(pdb_path, pdb_file), name="co_crystal")
protein, rna = splitEntities(structure) # split RNA and protein from structure

ss = getSS(prefix)
# print(ss)

with open("{}/{}-dssr.json".format(pdb_path, prefix)) as FH:
   data = json.load(FH, object_pairs_hook=collections.OrderedDict) 


protein_interactions,ss_dict = getInteractions(protein, rna, prefix)
pairs,backbone_edges, interaction_edges, interaction_types, stacks = getEdges(data, protein_interactions, ss_dict)
lw_values = getLW(data)
print(lw_values)

#update: added functions to extract all H-bond interactions from dssr and to add H-bond labels to interaction_types object
hbond_set = hbondExtractor(data)
interaction_types  = labelHbondEdges(interaction_types, hbond_set)

all_edges = pairs + backbone_edges + interaction_edges + stacks

d3 = d3graph(support=None, collision=0.5)
df = pd.DataFrame(all_edges, columns=['source', 'target'])
df['weight'] = [100]*len(pairs) + [100]*(len(backbone_edges)) + [5]*(len(interaction_edges)) + [20]*(len(stacks))
adjmat = vec2adjmat(df['source'], df['target'], weight=df['weight'])

###############################################################

d3.graph(adjmat)

d3.set_edge_properties(directed=True) # setting earlier to then update?

chains_list, centroid_rnaprodb_map, rotationMatrix, centroids_3d = getChainsAndPca(structure, interaction_edges)



d3.node_properties = processNodes(d3.node_properties)
ADD_PCA = True
if(ADD_PCA):
   d3.node_properties = addPcaToGraph(d3.node_properties, centroid_rnaprodb_map, centroids_3d)

##ADD RNAscape and ViennaRNA
try:
   d3.node_properties = addRNAscapeToGraph(d3.node_properties, structure, data, prefix)
except Exception as e:
   pass
try:
   d3.node_properties = addViennaToGraph(d3.node_properties, data, prefix)
except Exception as e:
   pass
#print("".join(["\n"]*100))

d3.edge_properties = processEdges(d3.edge_properties, backbone_edges, stacks, pairs, interaction_types, centroids_3d)

# d3.show(filepath='{}/output/{}.html'.format(home, pdb_file), show_slider=False, showfig=False)
# click={'fill': None, 'stroke': '#F0F0F0', 'size': 2.5, 'stroke-width': 10} # add inside d3 show to highlight click
final_json = d3.show(filepath='{}/output/{}.html'.format(home, pdb_file), show_slider=False, showfig=False)
# print(final_json)
final_json_object = json.loads(final_json)
ss_json = processSS(ss)
final_json_object["ss"] = ss_json
final_json_object["chainsList"] = chains_list
final_json_object["rotationMatrix"] = rotationMatrix.tolist() # used to orient NGLViewer camera to the PCA
final_json_object["tooLarge"] = TOO_LARGE

for edge in final_json_object["links"]:
    edge_tuple = (edge["source_id"], edge["target_id"])
    print(edge_tuple)
    print(lw_values)
    if edge_tuple in lw_values:
        print("FOUND " + str(edge_tuple))
        edge["LW"] = lw_values[edge_tuple]
    else:
        edge["LW"] = None

for edge in final_json_object["links"]:
   del edge['weight']
   del edge['weight_scaled']
   del edge['edge_distance']
   del edge['edge_weight']

final_json_object = json.dumps(final_json_object)
print(final_json_object)

with open('{}/output/{}_graph.json'.format(home, prefix), 'w') as outfile:
    outfile.write(final_json_object)
