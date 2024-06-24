from Bio import PDB
import numpy as np
import RNA
def getPairfromDBN(dbn="....((((....))..))...(())...&....((..)....)...()..(((..))).."):# seq="ABCDEFGHIJKLMNOPQRSTUVWXYZ12"):
    try:
        dbn = "".join(dbn)
    except:
        pass
    stack = []
    pairs = []
    dbn = dbn.replace("&","")
    for i in range(len(dbn)):
        if dbn[i] in [".","[","]"] :
            continue
        if dbn[i] == "(":
            stack.append(i)
        if dbn[i] == ")":
            if len(stack) == 0:
                continue
            pairs.append((stack[-1],i))
            stack = stack[:-1]
    return pairs, stack


def rot2eul(rotation_matrix):
    ### first transform the matrix to euler angles
    r =  Rotation.from_matrix(rotation_matrix)
    angles = r.as_euler("zyx",degrees=False)
    return angles

def addViennaToGraph(node_properties, edge_properties, data, prefix, scale = 5):
    dbn = []
    rnaprodb_ids = []
    for nt in data['nts']:
        rnaprodb_ids += ["{}:{}:{}".format(nt["chain_name"],nt["nt_resnum"],nt['nt_id'].split(".")[-1])]
        dbn += [nt["dbn"]]
    coords = RNA.get_xy_coordinates("".join(dbn))

    ### set SS pair keys
    ss_pair_idxs, _ = getPairfromDBN(dbn)
    ss_pairs = []
    for item in ss_pair_idxs:
        ss_pairs += [(rnaprodb_ids[item[0]], rnaprodb_ids[item[1]])]
        ss_pairs += [(rnaprodb_ids[item[1]], rnaprodb_ids[item[0]])]
    ####################
    
    centers = []
    for node_id, node in node_properties.items():
        if node['rnaprodb_id'] in rnaprodb_ids: # node has a centroid computed for it!
            idx = rnaprodb_ids.index(node['rnaprodb_id'])
            node['viennarna_x'] = coords.get(idx).X*scale
            node['viennarna_y'] = coords.get(idx).Y*scale
            #node['x'] = coords.get(idx).X*scale
            #node['y'] = coords.get(idx).Y*scale
            centers.append([node['viennarna_x'],node['viennarna_y']])
    
    centers = np.array(centers)
    centroid = np.mean(centers, axis=0)
    bounds = np.array([np.max(centers[:,0]), np.min(centers[:,0]), np.max(centers[:,1]),
        np.min(centers[:,1])])
    pro_updated = []
    shuf_edges = list(edge_properties.keys())
    np.random.seed(0)
    np.random.shuffle(shuf_edges)
    
    c = 0
    f = 0
    for edge_id in shuf_edges:
        c+=1
        edge = edge_properties[edge_id]
        edge['tertiary'] = False
        if edge['my_type'] in ["pair", "none"]:
            if edge['source_type'] == 'n' and edge['target_type'] == 'n':
                cand = (edge['source_id'], edge['target_id'])
                if cand not in ss_pairs:
                    f += 1
                    edge['tertiary'] = True


        if edge['source_id'] in rnaprodb_ids and edge['target_id'] not in rnaprodb_ids:
            idx = 0
        elif edge['target_id'] in rnaprodb_ids and edge['target_id'] not in rnaprodb_ids:
            idx = 1
        else:
            continue
        if edge_id[1-idx] in pro_updated:
            continue
        nuc_coords = np.array([node_properties[edge_id[idx]]['viennarna_x'],
            node_properties[edge_id[idx]]['viennarna_y']])
        try:
            prot_coords = np.array([node_properties[edge_id[1-idx]]['x'], node_properties[edge_id[1-idx]]['y']])
        except:
            continue
        
        def force_bound(coord, bounds, relax = 1.1):
            buonds = bounds*relax
            if coord[0] < bounds[0] and coord[0] > bounds [1] and coord[1] < bounds[2] and coord[1] > bounds [3]:
                return coord
            else:
                while not(coord[0] < bounds[0] and coord[0] > bounds [1] and coord[1] < bounds[2] and coord[1] > bounds [3]):
                    coord = coord * 0.9
            return coord
        if np.linalg.norm(nuc_coords - prot_coords) > 20:
            anchor = nuc_coords# or nuc_coords
            prot_coords = anchor + (20*(anchor - prot_coords)/np.linalg.norm(anchor -prot_coords))
        #prot_coords = force_bound(prot_coords, bounds)
        node_properties[edge_id[1-idx]]['viennarna_x'] = prot_coords[0]
        node_properties[edge_id[1-idx]]['viennarna_y'] = prot_coords[1]
        pro_updated.append(edge_id[1-idx])        
    #print(pro_updated, len(pro_updated))
    return node_properties

if __name__ == "__main__":
    print(getPairfromDBN())

