def hbondExtractor(data) -> list:
    all_pair_info = set()
    for bond in data['hbonds']:
        atom1 = bond["atom1_id"].split("..")[1]
        atom2 = bond["atom2_id"].split("..")[1]
        pair_info = (atom1, atom2)
        pair_info_rev = (atom2, atom1)
        if pair_info not in all_pair_info:
            all_pair_info.add(pair_info)
            all_pair_info.add(pair_info_rev)
    return all_pair_info

def labelHbondEdges(interaction_types_object, hbonds_set, ss_dict) -> dict:
    for hbond in hbonds_set:
        spl1 = hbond[0].split(".")
        spl2 = hbond[1].split(".")
    for pair, interactions in interaction_types_object.items():
        refined_pair = []
        for residue in pair:
            residue_type = residue.split(":")[1]
            residue_number = residue.split(":")[2] 
            refined_pair.append(f"{residue_type}.{residue_number}.")
        for item in hbonds_set:
            if refined_pair[0] in list(item)[0] and refined_pair[1] in  list(item)[1]:
                interactions.add('hbond')
    return interaction_types_object

