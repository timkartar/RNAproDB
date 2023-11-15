from utilities import node_to_text, parse_node, d3to1

nt_colors = {'A': '#FF9896',#'#90cc84',
    'C': '#DBDB8D',#'#AEC7E8',
    'G': '#90cc84',#'#DBDB8D',
    'U': '#AEC7E8'#'#FF9896'
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
        try:
            name = parsed_node[1]
        except:
            print(parsed_node, node)
        pos = str(parsed_node[2])
        chain = parsed_node[3]

        # global changes
        node_properties[node]['size'] = 25 # original 20
        node_properties[node]['opacity']= 1
        node_properties[node]['edge_size']= 1 # original 5
        node_properties[node]['fontcolor']= 'black'

        #ID of chain:number
        node_properties[node]['rnaprodb_id'] = "{}:{}".format(chain, pos)

        if(parsed_node[0] == 'n'): # is a nucleotide
            node_properties[node]['color']= nt_colors[name] #use nt color scheme
            tooltip = 'Nucleotide: ' + name +"\nPosition: " + pos + "\nChain: " + parsed_node[3]
            node_properties[node]['shape'] = 'circle' #is detected in  d3graphscript.js
            node_properties[node]['label']= name # empty label, use tooltip instead
            node_properties[node]['fontsize']= 25

        else: # is protein residue
            # Get one letter code
            if name in d3to1:
                one_letter_code = d3to1[name]
            else:
                one_letter_code = "X"
            node_properties[node]['fontsize']= 15
            ss = parsed_node[4]
            node_properties[node]['size'] = 30 # original 5
            node_properties[node]['color']= '#c6c6c6' #use gray by default
            node_properties[node]['label']= one_letter_code
            node_properties[node]['shape'] = 'rect' # detect in the JS
            if(ss == "H"): #Helix
                node_properties[node]['color']= 'white'
            elif(ss == "S"): #sheet
                node_properties[node]['color']= '#e8e8e8'
            elif(ss == "Unknown"): #Unknown
                node_properties[node]['color']= 'red'
            tooltip = 'Residue: ' + name +"\nPosition: " + pos + "\nChain: " + chain + "\nSec. Structure: " + ss
        node_properties[node]['tooltip']= tooltip
    return node_properties

def check_wc_pairing(edge_tuple):
    item1 = edge_tuple[0].split(":")[1]
    item2 = edge_tuple[1].split(":")[1]
    wc_pairs = ['AU','GC','UA','CG']
    if (item1 + item2) in wc_pairs:
        return True
    return False

def processEdges(edge_properties, backbone_edges, stacks, pairs, interaction_types):
    for edge in edge_properties:
        first_node,sec_node = parse_edge(edge)

        #global ids for each
        edge_properties[edge]['source_id'] = "{}:{}".format(first_node[3], first_node[2]) # chain:#
        edge_properties[edge]['target_id'] = "{}:{}".format(sec_node[3], sec_node[2]) # chain:#


        edge_properties[edge]['my_type'] = 'none'
        edge_properties[edge]['marker_end'] = ''
        edge_tuple = (node_to_text(first_node),node_to_text(sec_node)) #turn back into text to compare to backbone edge
        if edge_tuple in backbone_edges: #is a backbone edge. NOTE change to iterate through backbone edges instead!
            edge_properties[edge]['marker_start'] = ''
            # edge_properties[edge]['marker_end'] = 'arrow' # already set in set edge properties
            edge_properties[edge]['color'] = '#000000'#'#605f5f' # works!
            # edge_properties[edge]['label_color'] = 'red'
            # edge_properties[edge]['label_fontsize'] = 8
            # edge_properties[edge]['marker_color'] = 'red' # BROKEN!
            edge_properties[edge]['my_type'] = 'backbone'
            edge_properties[edge]['marker_end'] = 'arrow'
            edge_properties[edge]['edge_width'] = 6
        elif edge_tuple in stacks:
            edge_properties[edge]['marker_start'] = ''
            # edge_properties[edge]['my_type'] = 'link'
            # edge_properties[edge]['marker_end'] = 'square'
            edge_properties[edge]['edge_width'] = 2
        
        elif edge_tuple in pairs:
            if not check_wc_pairing(edge_tuple):
                edge_properties[edge]['color'] = '#F2936D'
            # edge_properties[edge]['marker_end'] = 'square'
            edge_properties[edge]['edge_width'] = 6
            edge_properties[edge]['my_type'] = 'pair'
            # edge_properties[edge]['my_type'] = "link"
        else:
            edge_properties[edge]['edge_width'] = 2
            edge_properties[edge]['marker_end'] = ''
            # edge_properties[edge]['my_type'] = "link"
            # edge_properties[edge]['directed'] = False
        
        ### segregate protein-RNA interaction edges
        if edge_tuple in interaction_types.keys():
            types = list(interaction_types[edge_tuple])
            if "major" in types and "minor" not in types:
                edge_properties[edge]['color'] = '#5DE3BA'
                edge_properties[edge]['edge_width'] = 4
            elif "minor" in types and "major" not in types:
                edge_properties[edge]['color'] = '#CE5DE3'
                edge_properties[edge]['edge_width'] = 4
            elif "major" in types and "minor" in types:
                edge_properties[edge]['color'] = 'black'
            elif "other" in types:
                edge_properties[edge]['color'] = 'black'
            else: #only backbone
                edge_properties[edge]['color'] = 'black' #for now
            if "hbond" in types:
                edge_properties[edge]['my_type'] = 'protein_rna_hbond'


    return edge_properties
