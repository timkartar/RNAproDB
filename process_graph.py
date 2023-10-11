from utilities import node_to_text, parse_node

nt_colors = {'A': '#00994C',
    'C': '#000099',
    'G': '#FBA922',
    'U': '#990000'
}

"""
Returns (('p'/'n', name, position, chain),('p'/'n', name, position, chain))
ASSUMES: proteins are represented by 3 letter code, and nucleotides by one letter!
"""
def parse_edge(node_id):
    first_node = parse_node(node_id[0])
    sec_node = parse_node(node_id[1])
    return first_node,sec_node

def processNodes(node_properties):
    for node in node_properties:
        parsed_node = parse_node(node) # ('p'/'n', name, position, chain, ss(protein only))
        name = parsed_node[1]
        pos = str(parsed_node[2])
        chain = parsed_node[3]
        
        # global changes
        node_properties[node]['size'] = 12 # original 20
        node_properties[node]['label']= "" # empty label, use tooltip instead
        node_properties[node]['opacity']= 0.705
        node_properties[node]['fontsize']= 10
        node_properties[node]['edge_size']= 1 # original 5

        if(parsed_node[0] == 'n'): # is a nucleotide
            node_properties[node]['color']= nt_colors[name] #use nt color scheme
            tooltip = 'Nucleotide: ' + name +"\nPosition: " + pos + "\nChain: " + parsed_node[3]
            node_properties[node]['fontcolor']= nt_colors[name]
        else: # is protein residue
            ss = parsed_node[4]
            node_properties[node]['size'] = 5 # original 5
            node_properties[node]['color']= 'gray' #use gray by default
            node_properties[node]['label']= name + pos # empty label, use tooltip instead
            node_properties[node]['fontcolor']= 'black'
            if(ss == "H"): #Helix
                node_properties[node]['color']= 'white'
            elif(ss == "S"): #Helix
                node_properties[node]['color']= 'black'
            elif(ss == "Unknown"): #Unknown
                node_properties[node]['color']= 'red'
            tooltip = 'Residue: ' + name +"\nPosition: " + pos + "\nChain: " + chain + "\nSec. Structure: " + ss
        node_properties[node]['tooltip']= tooltip
    return node_properties

def processEdges(edge_properties, backbone_edges, stacks):
    for edge in edge_properties:
        first_node,sec_node = parse_edge(edge)
        edge_tuple = (node_to_text(first_node),node_to_text(sec_node)) #turn back into text to compare to backbone edge
        if edge_tuple in backbone_edges: #is a backbone edge. NOTE change to iterate through backbone edges instead!
            edge_properties[edge]['marker_start'] = ''
            # edge_properties[edge]['marker_end'] = 'arrow' # already set in set edge properties
            edge_properties[edge]['color'] = 'red' # works!
            # edge_properties[edge]['label_color'] = 'red'
            # edge_properties[edge]['label_fontsize'] = 8
            # edge_properties[edge]['marker_color'] = 'red' # BROKEN!
        if edge_tuple in stacks:
            edge_properties[edge]['marker_start'] = ''
            # edge_properties[edge]['marker_end'] = 'square'
        else:
            edge_properties[edge]['marker_end'] = ''
            # edge_properties[edge]['directed'] = False
    return edge_properties