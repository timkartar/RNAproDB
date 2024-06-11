from Bio import PDB
import numpy as np
import RNA

def rot2eul(rotation_matrix):
    ### first transform the matrix to euler angles
    r =  Rotation.from_matrix(rotation_matrix)
    angles = r.as_euler("zyx",degrees=False)
    return angles

def addViennaToGraph(node_properties, data, prefix):
    dbn = []
    rnaprodb_ids = []
    for nt in data['nts']:
        rnaprodb_ids += ["{}:{}".format(nt["chain_name"],nt["nt_resnum"])]
        dbn += [nt["dbn"]]
    coords = RNA.get_xy_coordinates("".join(dbn))
    for node_id, node in node_properties.items():
        if node['rnaprodb_id'] in rnaprodb_ids: # node has a centroid computed for it!
            idx = rnaprodb_ids.index(node['rnaprodb_id'])
            node['viennarna_x'] = coords.get(idx).X
            node['viennarna_y'] = coords.get(idx).Y
    return node_properties

