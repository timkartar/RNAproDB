"""
Take selection of nodes as input
Process the JSON we already have to only include certain residues
Can nx process this, and easily find subgraphs?
Start with pre-computed JSON first, JSON dumps to file
"""

"""
Search for JSON in output folder
"""
import os, sys
import json
import networkx as nx

home =  os.path.dirname(os.path.abspath(__file__)) #change this line only 

pdb_path = "{}/output/".format(home)

# for use in generating source and target labels
def create_node_index_mapping(nodes):
    return {node: index for index, node in enumerate(nodes)}


def readJSON(pdbid):
    with open('{}/output/{}.json'.format(home, pdbid), 'r') as infile:
        data_string = infile.read()
    first_encode_data = json.loads(data_string)
    data = json.loads(first_encode_data)
    return data

def json_to_nx(data):
    # Create a directed graph using NetworkX
    G = nx.DiGraph()
    
    # Add nodes to the graph
    for node in data['nodes']:
        node_id = node['rnaprodb_id']
        del node['rnaprodb_id']
        G.add_node(node_id, **node)
    
    for link in data['links']:
        source = link['source_id']  # Using 'source_id' as identifier
        target = link['target_id']  # Using 'target_id' as identifier
        # print(f"Adding edge: {source} -> {target}")  # Debug line
        # Remove source and target ids to store all other attributes
        del link['source_id'], link['target_id']
        G.add_edge(source, target, **link)

    return G

def get_neighbors_within_distance(G, nodes, distance=2):
    nodes_set = set(nodes)
    for _ in range(distance):
        for node in list(nodes_set):
            nodes_set = nodes_set.union(set(G.neighbors(node)))
    return nodes_set

def nx_to_json(G, node_index_mapping):
    data = {
        'nodes': [{'rnaprodb_id': node, **G.nodes[node]} for node in G.nodes()],
        'links': [
            {
                'source': node_index_mapping[u],
                'target': node_index_mapping[v],
                **G[u][v]
            }
            for u, v in G.edges() if u in node_index_mapping and v in node_index_mapping
        ]
    }
    return data
def update_subgraph_edge_indices(subgraph, node_index_mapping):
    for u, v, edge_data in subgraph.edges(data=True):
        edge_data['source'] = node_index_mapping[u]
        edge_data['target'] = node_index_mapping[v]

if __name__ == "__main__":
    data = readJSON(sys.argv[1])
    G = json_to_nx(data)
    user_nodes = ["C:955", "A:646"]
    #user_nodes = ["C:6", "J:25"]
    subgraph_nodes = get_neighbors_within_distance(G,user_nodes)
    # print(subgraph_nodes)
    subgraph_edges = [(u, v) for u, v in G.edges() if u in subgraph_nodes and v in subgraph_nodes]
    # print(subgraph_edges) 
    subgraph = G.subgraph(subgraph_nodes).edge_subgraph(subgraph_edges)
    # print(subgraph.edges())
    node_index_mapping = create_node_index_mapping(subgraph.nodes())
    # print(subgraph['C:917']['C:918'])
    update_subgraph_edge_indices(subgraph,node_index_mapping)
    # print(subgraph['C:917']['C:918'])
    # print(node_index_mapping)
    subgraph_json = nx_to_json(subgraph, node_index_mapping)
    # print(subgraph_json)
    final_json_str = json.dumps(subgraph_json)
    final_json_str = json.dumps(final_json_str)
    print(final_json_str)
