from Bio import PDB
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation   

def rot2eul(rotation_matrix):
    ### first transform the matrix to euler angles
    r =  Rotation.from_matrix(rotation_matrix)
    angles = r.as_euler("zyx",degrees=False)
    return angles

def getChainsAndPca(structure, interaction_edges):
    chains_list = []
    all_centroids = []
    rna_centroids = []
    centroids_3d = {}

    aa_set = set() # which amino acids interact
    for edge in interaction_edges:
        split_edge = edge[1].split(":") #('C:C:973', 'A:LEU:278:H')
        aa_set.add("{}:{}".format(split_edge[0], split_edge[2]))


    for model in structure: # assume one model since bio assembly
        for chain in model:
            chain_dict = {} # wrapper holding resiudes and ID
            residue_list = [] # holds residues (name, pos, chain)
            
            chain_name = chain.get_id()
                        
            for residue in chain:
                residue_dict = {}
                residue_id = residue.get_id()[1]  # Gets the residue sequence number
                residue_name = residue.get_resname()

                residue_dict["name"] = residue_name
                residue_dict["pos"] = residue_id
                residue_dict["chain"] = chain_name
                residue_dict["is_aa"] = PDB.is_aa(residue)

                rnaprodbid = "{}:{}".format(residue_dict["chain"], residue_dict['pos'])


                if not residue_dict["is_aa"]:
                # if rnaprodbid in aa_set:
                    # Code below is to calculate centroid for use in PCA
                    atoms = [atom.get_vector() for atom in residue]
                    atom_coords = np.array([list(atom) for atom in atoms])
                    # Compute the mean (x,y,z) - the centroid of the residue
                    centroid = np.mean(atom_coords, axis=0)
                    rna_centroids.append(centroid)
                    all_centroids.append(centroid)
                    centroids_3d[rnaprodbid] = centroid
                elif rnaprodbid in aa_set:
                    atoms = [atom.get_vector() for atom in residue]
                    atom_coords = np.array([list(atom) for atom in atoms])
                    # Compute the mean (x,y,z) - the centroid of the residue
                    centroid = np.mean(atom_coords, axis=0)
                    all_centroids.append(centroid)
                    centroids_3d[rnaprodbid] = centroid

                residue_list.append(residue_dict)

            chain_dict["chainId"] = chain_name
            chain_dict["residues"] = residue_list
            chains_list.append(chain_dict)

    rna_centroids_array = np.array(rna_centroids)
    all_centroids_array = np.array(all_centroids)
    # Perform PCA to reduce to 2D
    SCALAR = 20
    
    pca = PCA(n_components=2)
    pca.fit(rna_centroids_array)
    reduced_centroids = pca.transform(np.array(all_centroids)) * SCALAR
    
    ## Alternative way that returns the rotation matrix
    ## We can send the rotation matrix to NGL viewer to set 
    ## the camera angle matching the explorer.
    ## var m = stage.viewerControls.getOrientation()
    ## stage.viewerControls.orient(m)
    '''
    U, S, Vt = np.linalg.svd(rna_centroids_array)
    v_rotated = all_centroids_array @ Vt.T 

    # Vt.T is the rotation matrix and its inverse is Vt (may
    # need either one)
    
    reduced_centroids = v_rotated[:,:2] *SCALAR # first two dimensions
    '''
    Vt = np.array([[0]*3]*3)
    # reduced_centroids[:,0] = -1*reduced_centroids[:,0]

    reduced_centroids = reduced_centroids - np.mean(reduced_centroids, axis=0)
    # print(reduced_centroids)
    # add the pca to chains_list
    i = 0
    # # Plot the reduced centroids
    # plt.scatter(reduced_centroids[:, 0], reduced_centroids[:, 1])
    # plt.xlabel('Principal Component 1')
    # plt.ylabel('Principal Component 2')
    # plt.title('PCA of Residue Centroids')
    # plt.show()
    centroid_rnaprodb_map = {}


   

    for chain in chains_list:
        for residue in chain["residues"]:
            rnaprodbid = "{}:{}".format(residue["chain"], residue['pos'])
            if not residue["is_aa"] or rnaprodbid in aa_set: # only nts and interacting residues
                residue["xpos"] = reduced_centroids[i][0]
                residue["ypos"] = reduced_centroids[i][1]
                centroid_rnaprodb_map[rnaprodbid] = reduced_centroids[i]
                i += 1
    # print("I IS ")
    # print(i)
    # print("Length of centroids")
    # print(len(reduced_centroids))
    # exit()
    r = rot2eul(Vt.T)
    r = np.array(Vt.T.tolist() + [r.tolist()])
    return chains_list, centroid_rnaprodb_map, r, centroids_3d
    # return chains_list, centroid_rnaprodb_map, Vt


def addPcaToGraph(node_properties, centroid_rnaprodb_map, centroids_3d):
    print(node_properties)
    for node_id, node in node_properties.items():
        print(node)
        print(node['rnaprodb_id'])
        if node['rnaprodb_id'] in centroid_rnaprodb_map: # node has a centroid computed for it!
            node['x'] = centroid_rnaprodb_map[node['rnaprodb_id']][0]
            node['y'] = centroid_rnaprodb_map[node['rnaprodb_id']][1]
            node['3dx'] = centroids_3d[node['rnaprodb_id']][0]
            node['3dy'] = centroids_3d[node['rnaprodb_id']][1]
            node['3dz'] = centroids_3d[node['rnaprodb_id']][2]
    return node_properties

