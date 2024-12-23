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


def readJSON(pdbid, algorithm):
    # with open('{}/output/{}_graph.json'.format(home, pdbid), 'r') as infile:
    #     data_string = infile.read()
    with open("{}/output/{}_{}_graph.json".format(home,pdbid, algorithm), 'r') as infile:
        json_output = json.load(infile)
        return json_output

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

def get_neighbors_within_distance(G, nodes, distance=1):
    nodes_set = set(nodes)
    for _ in range(distance):
        for node in list(nodes_set):
            nodes_set = nodes_set.union(set(G.neighbors(node))) # will break if node not in graph!!! Handle this later
            nodes_set = nodes_set.union(set(G.predecessors(node)))
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
    data = readJSON(sys.argv[1], sys.argv[3]) # assumes pre-computed JSON file, reads in pdbid
    G = json_to_nx(data) # convert json graph to nx
    user_nodes = sys.argv[2].split(",")
    user_nodes = list(filter(None, user_nodes)) #filter empty nodes
    

    subgraph_nodes = get_neighbors_within_distance(G,user_nodes) #list of nodes that should be in subgraph
    subgraph_edges = [(u, v) for u, v in G.edges() if u in subgraph_nodes and v in subgraph_nodes] # get the edges
    subgraph = G.subgraph(subgraph_nodes).edge_subgraph(subgraph_edges) #generate the new subgraph
    
    # each node for d3 must have a unique numerical index
    node_index_mapping = create_node_index_mapping(subgraph.nodes()) # compute the unique mapping for each subgraph node
    update_subgraph_edge_indices(subgraph,node_index_mapping) # update the graph to have a new source/target label

    subgraph_json = nx_to_json(subgraph, node_index_mapping) # convert nx back to JSON
    final_json_str = json.dumps(subgraph_json)
    print(final_json_str) #use console output
