import json
def processWH(whbond_data):
    output = {}
    for water in whbond_data:
        for nt_data in whbond_data[water][0]:
            for aa_data in whbond_data[water][1]:
                nt_part= nt_data[1]
                aa_part = aa_data[1]
                edge_id = (nt_part['rnaprodb_id'], aa_part['rnaprodb_id'])
                if edge_id not in output:
                    output[edge_id] = []
                toappend = {}
                toappend['water'] = water
                toappend['nt_role'] = nt_part['type']
                toappend['aa_role'] = aa_part['type']
                toappend['nt_distance'] = str(nt_part['distance_to_water']) + " Angstrom"
                toappend['aa_distance'] = str(aa_part['distance_to_water']) + " Angstrom"
                toappend['nt_atom'] = nt_part['atom']
                toappend['aa_atom'] = aa_part['atom']
                output[edge_id].append(toappend)

    # weiyu: add rna-rna water mediated hbonds
    # print(whbond_data['S-2:106:'])
    # print("after")

    for water in whbond_data:
        nt_data_list = whbond_data[water][0]

        # weiyu: need to include both directions: source --> target and target --> source
        for i in range(len(nt_data_list)):
            for j in range(len(nt_data_list)):
                if i == j:
                    continue
                nt1_part = nt_data_list[i][1]
                nt2_part = nt_data_list[j][1]
                edge_id = (nt1_part['rnaprodb_id'], nt2_part['rnaprodb_id'])

                if edge_id not in output:
                    output[edge_id] = []

                toappend = {}
                toappend['water'] = water
                toappend['nt_role'] = nt1_part['type']
                toappend['aa_role'] = nt2_part['type']
                toappend['nt_distance'] = str(nt1_part['distance_to_water']) + " Angstrom"
                toappend['aa_distance'] = str(nt2_part['distance_to_water']) + " Angstrom"
                toappend['nt_atom'] = nt1_part['atom']
                toappend['aa_atom'] = nt2_part['atom']

                output[edge_id].append(toappend)

    return output

def makeTooltip(json_obj, ntss_dict, whbond_data, hbond_set, hbond_data):
    node_idxs = {}
    node_id_idxs = {}
    edge_idxs = {}
    edge_dict = {}
    c = 0

    whbonds = processWH(whbond_data)

    for node in json_obj['nodes']:
        node_idxs[node['name']] = c
        node_id_idxs[node['rnaprodb_id']] = c
        #for key in node:
        #    print(key, node[key])
        c+=1

    c= 0
    for edge in json_obj['links']:
        #for item in edge:
        #    print(item, edge[item])
        #exit()
        try:
            source_name = edge['source_label']
            target_name = edge['target_label']
        except:
            #for item in edge:
            #    print(item, edge[item])
            pass
            #exit()

        if source_name not in edge_dict:
            edge_dict[source_name] = [c]
        else:
            edge_dict[source_name].append(c)

        if source_name not in edge_dict:
            edge_dict[source_name] = [c]
        else:
            edge_dict[source_name].append(c)

        table = {}
        table["Node1"] = edge['source_id']
        table["Node2"] = edge['target_id']
        table["Centroid distance"] = str(edge['distance_3d']) + " Angstrom"
        if edge['my_type'] == 'none':
            my_type = "other"
        else:
            my_type = edge['my_type']
        if edge['is_whbond']:
            table["Attributes"] = ", ".join([my_type, "water-mediated H bond"])
        else:
            table["Attributes"] = ", ".join([my_type])
        if edge['my_type'] == "pair":
            if edge['color'] == '#4169E1':
                table['pair_type'] = "Watson-Crick"
            else:
                table['pair_type'] = "Mismatch"
        if "LW" in edge:
            table["Leontis-Westhof class"] = edge['LW']

        source_idx = node_id_idxs[edge['source_id']]
        target_idx = node_id_idxs[edge['target_id']]
        source=json_obj['nodes'][source_idx]
        target=json_obj['nodes'][target_idx]

        cand1 = ".".join([source['name'].split(":")[0], source['name'].split(":")[1],
            source['name'].split(":")[2], source['icode']])
        cand2 = ".".join([target['name'].split(":")[0], target['name'].split(":")[1],
            target['name'].split(":")[2], target['icode']])
        hbond_cand = (cand1, cand2)
        if hbond_cand in hbond_set:
            idx = 1
            for item in hbond_data[hbond_cand]:
                table["Hbond {}".format(idx)] = ""
                #print(item)
                for key, val in item.items():
                    table["\\t{}_{}".format(key, idx)] = str(val)
                idx+=1

        whbond_cand = (edge['source_id'],edge['target_id'])

        if whbond_cand in whbonds:
            idx = 1
            for item in whbonds[whbond_cand]:
                table["Water-mediated Hbond {}".format(idx)] = ""
                for key, val in item.items():
                    table["\\t{}_{}".format(key, idx)] = str(val)
                idx+=1
            #break
        edge['tooltip_table'] = json.dumps(table)
        c+=1
    
    for node in json_obj['nodes']:
        table = {}
        table['Name'] = node['name'].split(":")[1]
        table['Chain'] = node['name'].split(":")[0]
        if "Nucleotide" in node['node_tooltip']:
            table['Residue number'] = "".join(node['name'].split(":")[2:])
            try:
                table['Structural motif'] = ntss_dict[node['name']]
            except:
                pass
        else:
            table['Residue number'] = node['name'].split(":")[2]
            table['Structural motif'] = node['name'].split(":")[3]

        node['tooltip_table'] = json.dumps(table)
    
    #for edge in json_obj['links']:
    #    edge['tooltip_table'] = json.dumps({"wqwrqwrq":"qFFQ"})
        
    return json_obj
