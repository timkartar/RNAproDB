import json
import copy
import itertools
import numpy as np

def makeSSgraph(json_obj, dssrss):
    ret = {"nodes": [],
            "links": []
    }
    
    ssnodes = {}
    for key, val in dssrss.items():
        for item in val:
            itemstr = ":".join([str(i) for i in item])
            if itemstr not in ssnodes:
                ssnodes[itemstr] = [key]
            else:
                ssnodes[itemstr].append(key)

    ######## check ssnode connections #########
    json_edges = []
    for item in json_obj['links']:
        if item['my_type'] == 'backbone':
            json_edges.append((item['source_label'],item['target_label']))
    
    json_int_edges = []
    for item in json_obj['links']:
        if item['my_type'] not in[ 'backbone', 'pair']:
            json_edges.append((item['source_label'],item['target_label']))

    pos_dict = {}
    id_dict = {}
    for item in json_obj['nodes']:
        pos_dict[item['name']] = [float(item['viennarna_x']),float(item['viennarna_y'])]
        id_dict[item['name']] = item['rnaprodb_id']
 
    ssnodes_list = list(ssnodes.keys())
    for ssnode in ssnodes_list:
        coords = []
        members = []
        for item in ssnodes[ssnode]:
            try:
                coords.append(pos_dict[item])
                members.append(id_dict[item])
            except:
                continue
        if len(coords) == 0: 
            del ssnodes[ssnode]

        else:
            coords = np.mean(coords, axis=0)
            ret["nodes"].append({"id": ssnode, 'x':str(coords[0]), 'y':str(coords[1]),
                'members':members})
    
    opposite_nodes = {}
    for key, val in ssnodes.items():
        for item in val:
            if item not in opposite_nodes:
                opposite_nodes[item] = [key]
            else:
                opposite_nodes[item].append(key)

    ssedges = set()
    ssnodekeys = list(ssnodes.keys())
    for node1 in ssnodekeys:
        for node2 in ssnodekeys:
            if node1 == node2:
                continue
            cands  = list(itertools.product(ssnodes[node1], ssnodes[node2]))
            for cand in cands:
                if cand in json_edges:
                    ssedges.add((node1, node2))
                    ssedges.add((node2, node1))
                    break
    

    for item in json_obj['nodes']:
        if "Residue" in item['node_tooltip']:
            name = item['name']
            #print(name)
            #print(name, opposite_nodes)
            #exit()
            #if name not in opposite_nodes:
            #    continue
            hasEdge = False
            for ssnode in ssnodes:
                for json_node in ssnodes[ssnode]:
                    #print(ssnode, json_node)
                    cands  = [(name, json_node), (json_node, name)]
                    for cand in cands:
                        if cand in json_edges:
                            ssedges.add(("protein:" + name, ssnode))
                            ssedges.add((ssnode, "protein:" + name))
                            hasEdge = True
                        else:
                            pass
            if hasEdge:
                ret['nodes'].append({"id": "protein:" + name,
                    'x': item['viennarna_x'],
                    'y': item['viennarna_y']
                    })
        
    for ssedge in ssedges:
        ret["links"].append({"source": ssedge[0], "target": ssedge[1]})
    
    '''
    for item in json_obj['nodes']:
        name = item['name']
        if name in dssrss.keys():
    '''     
    
    
    return ret, opposite_nodes
