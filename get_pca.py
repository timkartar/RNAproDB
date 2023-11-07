from Bio import PDB
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

def getChainsAndPca(structure, interaction_edges):
    chains_list = []
    all_centroids = []

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


                if not residue_dict["is_aa"] or rnaprodbid in aa_set:
                # if rnaprodbid in aa_set:
                    # Code below is to calculate centroid for use in PCA
                    atoms = [atom.get_vector() for atom in residue]
                    atom_coords = np.array([list(atom) for atom in atoms])
                    # Compute the mean (x,y,z) - the centroid of the residue
                    centroid = np.mean(atom_coords, axis=0)
                    all_centroids.append(centroid)
                

                residue_list.append(residue_dict)

            chain_dict["chainId"] = chain_name
            chain_dict["residues"] = residue_list
            chains_list.append(chain_dict)

    centroids_array = np.array(all_centroids)
    # Perform PCA to reduce to 2D
    pca = PCA(n_components=2)
    SCALAR = 20
    reduced_centroids = pca.fit_transform(centroids_array) * SCALAR

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
    return chains_list, centroid_rnaprodb_map

def addPcaToGraph(node_properties, centroid_rnaprodb_map):
    print(node_properties)
    for node_id, node in node_properties.items():
        print(node)
        print(node['rnaprodb_id'])
        if node['rnaprodb_id'] in centroid_rnaprodb_map: # node has a centroid computed for it!
            node['x'] = centroid_rnaprodb_map[node['rnaprodb_id']][0]
            node['y'] = centroid_rnaprodb_map[node['rnaprodb_id']][1]
    return node_properties

