from Bio import PDB
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation   
from dnaprodb_rnascape.rnascape import rnascape


def rot2eul(rotation_matrix):
    ### first transform the matrix to euler angles
    r =  Rotation.from_matrix(rotation_matrix)
    angles = r.as_euler("zyx",degrees=False)
    return angles

def addRNAscapeToGraph(node_properties, edge_properties, structure, data, prefix, scale=20):
    model = structure[0]
    rnascape_coords, markers, ids, chids, dssrids, dssrout, prefix = rnascape(prefix, model,
                    data, cond_bulging=False, mFIG_PATH="./rnascape/output/processed_images/", mDSSR_PATH =
                'x3dna-dssr')
    rnaprodb_ids = ["{}:{}".format(chids[i],ids[i][1]) for i in range(len(ids))]
    for node_id, node in node_properties.items():
        if node['rnaprodb_id'] in rnaprodb_ids: # node has a centroid computed for it!
            idx = rnaprodb_ids.index(node['rnaprodb_id'])
            node['rnascape_x'] = rnascape_coords[idx][0]
            node['rnascape_y'] = rnascape_coords[idx][1]
            #node['x'] = rnascape_coords[idx][0]*scale
            #node['y'] = rnascape_coords[idx][1]*scale
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
        #node_properties[edge_id[1-idx]]['x'] = prot_coords[0]
        #node_properties[edge_id[1-idx]]['y'] = prot_coords[1]
        node_properties[edge_id[1-idx]]['rnascape_x'] = prot_coords[0]
        node_properties[edge_id[1-idx]]['rnascape_y'] = prot_coords[1]
        
        print(node_properties[edge_id[1-idx]])

    return node_properties

