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

def addRNAscapeToGraph(node_properties, structure, data, prefix):
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
            # node['x'] = rnascape_coords[idx][0]
            # node['y'] = rnascape_coords[idx][1]
    return node_properties

