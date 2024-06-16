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
from get_whbonds import getWHbonds #runHBplus

parser = MMCIFParser(QUIET=True)
home =  os.path.dirname(os.path.abspath(__file__))
backend =  os.path.dirname(os.path.abspath(__file__))
frontend = backend + "/../rnaprodb_frontend/"

pdb_path = frontend + "public/cifs/".format(home)
# pdb_file = "8fvi-assembly1.cif"

if len(sys.argv) > 1:
   prefix = sys.argv[1]

pdb_file = "{}-assembly1.cif".format(prefix)

TOO_LARGE = False

#structure = StructureData(os.path.join(pdb_path, pdb_file), name="co_crystal")

structure = parser.get_structure(prefix, os.path.join(pdb_path, pdb_file))
protein, rna = splitEntities(structure) # split RNA and protein from structure
#print(water_hbonds)
#exit()

data = runDSSR(structure, quiet=True, prefix=prefix, tmpdir="")
#for residue in structure.get_residues():
#    print(residue.get_id())
ss = getSS(prefix, data)

#with open("{}/{}-dssr.json".format(pdb_path, prefix)) as FH:
#   data = json.load(FH, object_pairs_hook=collections.OrderedDict) 


protein_interactions,ss_dict = getInteractions(protein, rna, prefix)
#for item in protein_interactions:
#    #if "HIS" in protein_interactions[item]:
#    print(item, protein_interactions[item])
#exit()
pairs,backbone_edges, interaction_edges, interaction_types, stacks = getEdges(data, protein_interactions, ss_dict)
#for item in pairs:
#    print(item)
#    #if "PSU" in item[0] or "PSU" in item[1]:
#    #    print("WHOAAAA")
#exit()

#update: added functions to extract all H-bond interactions from dssr and to add H-bond labels to interaction_types object
hbond_set = hbondExtractor(data)
interaction_types  = labelHbondEdges(interaction_types, hbond_set)
water_hbonds, interaction_types, whbond_data = getWHbonds(pdb_path, "{}-assembly1".format(prefix), structure, ss_dict,
        interaction_types)

#print(water_hbonds)
all_edges = pairs + backbone_edges + interaction_edges + stacks + water_hbonds

#print(len(all_edges), len(list(set(all_edges)))) not same ??
#exit()
#for edge in all_edges:
#    print(edge)
#exit()
d3 = d3graph(support=None, collision=0.5)
df = pd.DataFrame(all_edges, columns=['source', 'target'])
df['weight'] = [100]*len(pairs) + [100]*(len(backbone_edges)) + [5]*(len(interaction_edges)) +[20]*(len(stacks)) + [20]*len(water_hbonds)
adjmat = vec2adjmat(df['source'], df['target'], weight=df['weight'])

#for item in all_edges:
#    print(item)
#    #if "PSU" in item[0] or "PSU" in item[1]:
#    #    print("WHOAAAA")
#exit()
###############################################################

d3.graph(adjmat)

d3.set_edge_properties(directed=True) # setting earlier to then update?

chains_list, centroid_rnaprodb_map, rotationMatrix, centroids_3d = getChainsAndPca(structure,
        interaction_edges + water_hbonds) # include only water mediated AAs in PCA



d3.node_properties = processNodes(d3.node_properties)
#for item in d3.node_properties:
#    print(item, d3.node_properties[item])
#exit()

ADD_PCA = True
if(ADD_PCA):
    d3.node_properties = addPcaToGraph(d3.node_properties, centroid_rnaprodb_map, centroids_3d)
d3.edge_properties = processEdges(d3.edge_properties, backbone_edges, stacks, pairs, interaction_types, centroids_3d)


edges = list(d3.edge_properties.keys())
for edge in edges:
    if edge[0] not in d3.node_properties.keys() or edge[1] not in d3.node_properties.keys() :
        del d3.edge_properties[edge]

nodes = list(d3.node_properties.keys())
for node in nodes:
    if "Residue" not in d3.node_properties[node]['tooltip']:
        continue
    has_edge = False
    for edge in d3.edge_properties:
        if node in edge:
            has_edge = True
            break
    if has_edge == False:
        del d3.node_properties[node]
##ADD RNAscape and ViennaRNA
d3.node_properties = addRNAscapeToGraph(d3.node_properties, d3.edge_properties, structure, data, prefix)

#for item in d3.node_properties:
#    try:
#        print(item, d3.node_properties[item]['rnascape_x'])
#    except:
#        print(item)
#exit()
d3.node_properties = addViennaToGraph(d3.node_properties, d3.edge_properties, data, prefix)

#coord_type = "viennarna" ## DUMMY REPLACE FOR TESTING / COMMENT OUT AND MAKE OPTION IN FRONTEND
#if coord_type == "pca":
for node in d3.node_properties:
    print(node, d3.node_properties[node])
    if 'x' not in d3.node_properties[node].keys():
        continue
    d3.node_properties[node]['pca_x'] = d3.node_properties[node]['x']
    d3.node_properties[node]['pca_y'] = d3.node_properties[node]['y']
    #for node in d3.node_properties:
    #    print("pca", d3.node_properties[node])
    #pass
'''
if coord_type == "rnascape":
    for node in d3.node_properties:
        d3.node_properties[node]['x'] = d3.node_properties[node]['rnascape_x']
        d3.node_properties[node]['y'] = d3.node_properties[node]['rnascape_y']
        #print("rnascape", d3.node_properties[node])
        #print(d3.node_properties[node])
#exit()
if coord_type == "viennarna":
    for node in d3.node_properties:
        d3.node_properties[node]['x'] = d3.node_properties[node]['viennarna_x']
        d3.node_properties[node]['y'] = d3.node_properties[node]['viennarna_y']

   
for node in d3.node_properties:
    print(node, d3.node_properties[node])

for edge in d3.edge_properties:
    print(edge, d3.edge_properties[edge])
'''

###################################################
#
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
final_json_object["whbond_data"] = whbond_data
# extract LW annotations
lw_values = getLW(data)

def flip_lw_value(lw):
    if len(lw) == 3:
        return lw[0] + lw[2] + lw[1]
    return lw

processed_edges = set()

for edge in final_json_object["links"]:
    if edge["my_type"] == "pair":
        source_num = int(edge["source_id"].split(':')[1])
        target_num = int(edge["target_id"].split(':')[1])

        if source_num < target_num:
            edge_tuple = (edge["source_id"], edge["target_id"])
            reverse_lw = False
        else:
            edge_tuple = (edge["target_id"], edge["source_id"])
            reverse_lw = True

        if edge_tuple in processed_edges:
            # print("Skipping already processed edge:", edge_tuple)
            continue

        processed_edges.add(edge_tuple)

        if edge_tuple in lw_values:
            lw_value = lw_values[edge_tuple]
            if reverse_lw:
                lw_value = flip_lw_value(lw_value)
            edge["LW"] = lw_value
      #       print("ADDING LW FOR:", edge_tuple, "WITH VALUE:", lw_value)
      #   else:
      #       print("LW NOT FOUND FOR", edge_tuple, "IN ORIGINAL ORDER")



for edge in final_json_object["links"]:
   del edge['weight']
   del edge['weight_scaled']
   del edge['edge_distance']
   del edge['edge_weight']

final_json_object = json.dumps(final_json_object)
print(final_json_object)

with open('{}/output/{}_graph.json'.format(home, prefix), 'w') as outfile:
    outfile.write(final_json_object)
