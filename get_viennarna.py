from Bio import PDB
import numpy as np
import RNA

def rot2eul(rotation_matrix):
    ### first transform the matrix to euler angles
    r =  Rotation.from_matrix(rotation_matrix)
    angles = r.as_euler("zyx",degrees=False)
    return angles

def addViennaToGraph(node_properties, edge_properties, data, prefix, scale = 5):
    dbn = []
    rnaprodb_ids = []
    for nt in data['nts']:
        rnaprodb_ids += ["{}:{}".format(nt["chain_name"],nt["nt_resnum"])]
        dbn += [nt["dbn"]]
    coords = RNA.get_xy_coordinates("".join(dbn))
    for node_id, node in node_properties.items():
        if node['rnaprodb_id'] in rnaprodb_ids: # node has a centroid computed for it!
            idx = rnaprodb_ids.index(node['rnaprodb_id'])
            node['viennarna_x'] = coords.get(idx).X*scale
            node['viennarna_y'] = coords.get(idx).Y*scale
            #node['x'] = coords.get(idx).X*scale
            #node['y'] = coords.get(idx).Y*scale
    for edge_id, edge in edge_properties.items():
        print(edge_id, edge)
        if edge['source_id'] in rnaprodb_ids and edge['target_id'] not in rnaprodb_ids:
            idx = 0
        elif edge['target_id'] in rnaprodb_ids and edge['target_id'] not in rnaprodb_ids:
            idx = 1
        else:
            continue
        print(node_properties[edge_id[1-idx]])
        nuc_coords = np.array([node_properties[edge_id[idx]]['x'],
            node_properties[edge_id[idx]]['y']])
        prot_coords = np.array([node_properties[edge_id[1-idx]]['x'], node_properties[edge_id[1-idx]]['y']])

        if np.linalg.norm(nuc_coords - prot_coords) > 10*scale:
            prot_coords = nuc_coords + (10*(nuc_coords - prot_coords)/np.linalg.norm(nuc_coords -
                    prot_coords))
        node_properties[edge_id[1-idx]]['viennarna_x'] = prot_coords[0]
        node_properties[edge_id[1-idx]]['viennarna_y'] = prot_coords[1]
    
    return node_properties

