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
from get_interactions import getInteractions, getProteinSecStructure
from run_dssr import runDSSR
from get_edges import getEdges
from process_graph import processEdges, processNodes

parser = MMCIFParser()

#update: no need to change anymore
home =  os.path.dirname(os.path.abspath(__file__)) #change this line only 

pdb_path = "{}/pdb/".format(home)
# pdb_file = "8fvi-assembly1.cif"
prefix = '1ivs'
pdb_file = "{}-assembly1.cif".format(prefix)
original_structure = StructureData(os.path.join(pdb_path, pdb_file), name="co_crystal")
protein, rna = splitEntities(original_structure) # split RNA and protein from structure

# NOTE: here we save the cleaned temporary file before processing it in runDSSR. For larger files, this poses a challenge.
rna = cleanRNA(rna)
data = runDSSR(rna, quiet=True, prefix='1ivs')

protein_interactions = getInteractions(protein, rna, prefix)
pairs,backbone_edges, interaction_edges,stacks = getEdges(data, protein_interactions)
all_edges = pairs + backbone_edges + interaction_edges + stacks

df = pd.DataFrame(all_edges, columns=['source', 'target'])
df['weight'] = [20]*len(pairs) + [100]*(len(backbone_edges)) + [5]*(len(interaction_edges)) + [20]*(len(stacks))
d3 = D3Blocks()
d3.d3graph(df, filepath=None, showfig=False)
d3.D3graph.set_edge_properties(directed=True, marker_color='red') # setting earlier to then update?

# can probably pre-compute, then add to the dataframe and use that?
# iterate through nodes to change colors, label, etc.
d3.D3graph.node_properties = processNodes(d3.D3graph.node_properties)
d3.D3graph.edge_properties = processEdges(d3.D3graph.edge_properties, backbone_edges)

d3.D3graph.show(filepath='{}/output/{}.html'.format(home, pdb_file), click={'fill': None, 'stroke': '#F0F0F0', 'size': 2.5, 'stroke-width': 10})

# Logs the node data to the console when you click, in color_on_click() function!
# console.log("Clicked Node Data:", d3.select(this).data()[0]);
