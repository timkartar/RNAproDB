import json
import os
import numpy

def makeTable(json_output, path=None):
    data = dict(json_output)

    table = dict()
    table['Interaction pairs'] = ["Node 1, Node 2, Centroid distance, Attribute"]
    table['Hydrogen bonds'] = ["Node 1,Node 2,Distance,Atom1,Atom2"]
    table['Water-mediated hydrogen bonds'] = ["Water,Node 1,Node 2,NT_distance,AA_distance,NT_role,AA_role,NT_Atom,AA_Atom"]
    table['Base pairing info'] = ["Node 1,Node 2,Centroid distance, pair_type, Leontis-Westhof class"]
    table['Structural motif'] = ["Name,Chain,Number,Motif (DSSP annotation for protein/DSSR for RNA)"]

    for item in data['nodes']:
        tooltip = json.loads(item['tooltip_table'].replace("\\\\t",""))
        node1 = tooltip["Name"]
        chain = tooltip["Chain"]
        res = tooltip["Residue number"]
        try:
            struct = tooltip["Structural motif"]
        except:
            struct = "-"
        toappend = "{},{},{},{}".format(node1, chain, res, str(struct).replace("[","").replace("]","").replace("'","").replace(",",";"))
        table['Structural motif'].append(toappend)
    for item in data['links']:
        tooltip = json.loads(item['tooltip_table'].replace("\\\\t",""))
        node1 = tooltip["Node1"]
        node2 = tooltip["Node2"]
        distance = tooltip["Centroid distance"]
        attribute = tooltip["Attributes"]
        toappend = "{},{},{},{}".format(node1, node2, distance, attribute.replace(",",";"))
        table['Interaction pairs'].append(toappend)
        print(toappend)

        for key in tooltip:
            if "Hbond" in key and "Water" not in key:
                node1 = tooltip["Node1"]
                node2 = tooltip["Node2"]
                idx = key.split(" ")[1]
                distance = tooltip["distance_"+idx]
                nt_atom = tooltip["atom1_"+idx]
                aa_atom = tooltip["atom2_"+idx]
                toappend = "{},{},{},{},{}".format(node1, node2, distance, nt_atom, aa_atom)
                table['Hydrogen bonds'].append(toappend)
            elif  "Water" in key:
                node1 = tooltip["Node1"]
                node2 = tooltip["Node2"]
                idx = key.split(" ")[-1]
                water = tooltip["water_"+idx]
                nt_distance = tooltip["nt_distance_"+idx]
                aa_distance = tooltip["aa_distance_"+idx]
                nt_role = tooltip["nt_role_"+idx]
                aa_role = tooltip["aa_role_"+idx]
                nt_atom = tooltip["nt_atom_"+idx]
                aa_atom = tooltip["aa_atom_"+idx]
                toappend = "{},{},{},{},{},{},{},{},{}".format(water, node1, node2, nt_distance, aa_distance, nt_role, aa_role, nt_atom, aa_atom)
                table['Water-mediated hydrogen bonds'].append(toappend)
            elif "Leontis" in key:
                node1 = tooltip["Node1"]
                node2 = tooltip["Node2"]
                distance = tooltip["Centroid distance"]
                pair_type = tooltip["pair_type"]
                lw_class = tooltip["Leontis-Westhof class"]
                toappend = "{},{},{},{},{}".format(node1, node2, distance, pair_type, lw_class)
                table['Base pairing info'].append(toappend)
            
    #print(table['Water-mediated hydrogen bonds'])
    return table

if __name__ == "__main__":
    table  = makeTable("./output/1ivs_pca_graph.json")
    json.dump(table, open("1ivs_table.json","w"))
