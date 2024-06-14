from Bio import PDB
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation   
from dnaprodb_rnascape.rnascape import rnascape


def addRNAscapeToGraph(node_properties, edge_properties, structure, data, prefix, scale=30):
    np.random.seed(0)
    model = structure[0]
    rnascape_coords, markers, ids, chids, dssrids, dssrout, prefix = rnascape(prefix, model,
                    data, cond_bulging=False, mFIG_PATH="./rnascape/output/processed_images/", mDSSR_PATH =
                'x3dna-dssr')
    rnaprodb_ids = ["{}:{}".format(chids[i],ids[i][1]) for i in range(len(ids))]
    centers = []
    for node_id, node in node_properties.items():
        if node['rnaprodb_id'] in rnaprodb_ids: # node has a centroid computed for it!
            idx = rnaprodb_ids.index(node['rnaprodb_id'])
            node['rnascape_x'] = rnascape_coords[idx][0]*scale
            node['rnascape_y'] = rnascape_coords[idx][1]*scale
            #node['x'] = rnascape_coords[idx][0]*scale
            #node['y'] = rnascape_coords[idx][1]*scale
            if [node['rnascape_x'],node['rnascape_y']] in centers:
                rand = [np.random.random()*10,np.random.random()*10]
                node['rnascape_x'] += rand[0]
                node['rnascape_y'] += rand[1]
            centers.append([node['rnascape_x'],node['rnascape_y']])
    centers = np.array(centers)
    centroid = np.mean(centers, axis=0)
    bounds = np.array([np.max(centers[:,0]), np.min(centers[:,0]), np.max(centers[:,1]),
        np.min(centers[:,1])])

    pro_updated = []
    shuf_edges = list(edge_properties.keys())
    np.random.shuffle(shuf_edges)
    for edge_id in shuf_edges:
        edge = edge_properties[edge_id]
        #print(edge_id, edge)
        if edge['source_id'] in rnaprodb_ids and edge['target_id'] not in rnaprodb_ids:
            idx = 0
        elif edge['target_id'] in rnaprodb_ids and edge['target_id'] not in rnaprodb_ids:
            idx = 1
        else:
            continue
        #print(node_properties[edge_id[1-idx]])
        if edge_id[1-idx] in pro_updated:
            continue
        nuc_coords = np.array([node_properties[edge_id[idx]]['rnascape_x'],
            node_properties[edge_id[idx]]['rnascape_y']])
        prot_coords = np.array([node_properties[edge_id[1-idx]]['x'], node_properties[edge_id[1-idx]]['y']])
        
        def force_bound(coord, bounds, relax = 1.1):
            buonds = bounds*relax
            if coord[0] < bounds[0] and coord[0] > bounds [1] and coord[1] < bounds[2] and coord[1] > bounds [3]:
                return coord
            else:
                while not(coord[0] < bounds[0] and coord[0] > bounds [1] and coord[1] < bounds[2] and coord[1] > bounds [3]):
                    coord = coord * 0.9
            return coord

        if np.linalg.norm(nuc_coords - prot_coords) > 30:
            anchor = nuc_coords# or nuc_coords
            prot_coords = anchor + (30*(anchor - prot_coords)/np.linalg.norm(anchor -prot_coords))
        #prot_coords = force_bound(prot_coords, bounds)
        #node_properties[edge_id[1-idx]]['x'] = prot_coords[0]
        #node_properties[edge_id[1-idx]]['y'] = prot_coords[1]
        node_properties[edge_id[1-idx]]['rnascape_x'] = prot_coords[0]
        node_properties[edge_id[1-idx]]['rnascape_y'] = prot_coords[1]
        pro_updated.append(edge_id[1-idx])
        

    return node_properties

