# """
# IP: Returns if two nodes (tuple representation) are a backbone
# """
# def is_backbone_edge(node_1, node_2):
#     #if (both nucleotides) AND (one nt away from each other) AND (on same strand)
#     if(node_1[0] == 'nt' and node_2[0] == 'nt' and ((int(node_2[2]) - int(node_1[2]) == 1) or (int(node_2[2]) - int(node_1[2]) == -1)) and node_1[3] == node_2[3]):
#        return True
#     else:
#         return False
    
"""
Returns rnaprodb nucleotide text string (A:B:C) from dssr id representation
"""
def dssr_id_to_text(dssr_id):
    return ":".join(dssr_id.split(".")[2:-1])



"""
stack[nts_long]: '..C.G.901.,..C.A.972.'
"""
def get_stacking_interactions(dssr):
    stack_interactions = []
    for stack in dssr['stacks']:
        nts_long_split = stack['nts_long'].split(',')
        first_nucleotide = dssr_id_to_text(nts_long_split[0])
        sec_nucleotide = dssr_id_to_text(nts_long_split[1])
        stack_interactions.append((first_nucleotide, sec_nucleotide))
    return stack_interactions

"""
Returns base and backbone pairings from pre-processed (i.e., RNA only) DSSR 
"""
def getEdges(dssr, protein_interactions):
    pairs_dict = {}
    pairs=[]

    #add pairs to dictionary for easy lookup
    for pair in dssr["pairs"]:
        p1 = dssr_id_to_text(pair.get("nt1"))
        p2 = dssr_id_to_text(pair.get("nt2"))
        pairs_dict[p1] = p2
        pairs_dict[p2] = p1
    
    # add base pairing and self edges
    for nt in dssr["nts"]:
        nt_id = dssr_id_to_text(nt.get("nt_id"))
        if nt_id in pairs_dict:
            # base pairing edge
            pair_edge = (nt_id, pairs_dict[nt_id])
        else:
            # self edge
            pair_edge = (nt_id, nt_id)
        pairs.append(pair_edge)
    # add backbone edges
    backbone_edges = []
    for i, nt in enumerate(dssr["nts"]):
        nt_id = dssr_id_to_text(nt.get("nt_id"))
        try:
            nt_next = dssr_id_to_text(dssr["nts"][i+1].get("nt_id"))
        except Exception as e:
            print(nt_id, e)
            continue
        c1 = nt_id.split(":")[0]
        c2 = nt_next.split(":")[0]

        if c1 == c2:
            backbone_edges.append((nt_id, nt_next))
        
    interaction_edges = []
    for key, val in protein_interactions.items():
        for v in val:
            nt = ":".join((key.split(":")[:-1])) #be careful here
            ss = key.split(":")[3]
            interaction_edges.append((nt, v))
            
    interaction_edges = list(set(interaction_edges))
    # print(interaction_edges) #('C:C', 'A:PRO:277:H')
    stacks = get_stacking_interactions(dssr)

    return pairs,backbone_edges, interaction_edges,stacks