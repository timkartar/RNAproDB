from utilities import dssr_id_to_text, is_a_protein


"""
NOTE: There can be multiple stacking interactions, and stacks can include protein residues!
NOTE: Assuming stack is in order of the file, otherwise will need to figure out the order!
OrderedDict([('index', 12), ('num_nts', 7), ('nts_short', 'GAGAGAC'), ('nts_long', '1..C.G.942.,1..C.A.943.,1..C.G.944.,1..C.A.909.,1..C.G.945.,1..C.A.920.,1..C.C.947.')])
"""
def get_stacking_interactions(dssr, ss_dict):
    stack_interactions = []
    for stack in dssr['stacks']:
        nts_long_split = stack['nts_long'].split(',')
        num_nts = stack["num_nts"]
        for i,nt in enumerate(nts_long_split):
            if i != (num_nts - 1): # not last nt in stack, proceed
                first_nucleotide = dssr_id_to_text(nts_long_split[i])
                sec_nucleotide = dssr_id_to_text(nts_long_split[i+1])
                
                if(is_a_protein(first_nucleotide)):
                    protein_num = first_nucleotide.split(":")[2] # gets residue #!
                    ss = ss_dict[int(protein_num)]
                    first_nucleotide = first_nucleotide + ":{}".format(ss) # append ss to the protein residue
                if(is_a_protein(sec_nucleotide)):
                    protein_num = sec_nucleotide.split(":")[2] # gets residue #!
                    ss = ss_dict[int(protein_num)]
                    sec_nucleotide = sec_nucleotide + ":{}".format(ss) # append ss to the protein residue

                stack_interactions.append((first_nucleotide, sec_nucleotide))
    return stack_interactions

"""
Returns base and backbone pairings from pre-processed (i.e., RNA only) DSSR 
"""
def getEdges(dssr, protein_interactions, ss_dict):
    pairs_dict = {}
    pairs=[]
    backbone_edges = []
    num_nts = dssr["num_nts"]

    #add pairs to dictionary for easy lookup
    for pair in dssr["pairs"]:
        p1 = dssr_id_to_text(pair.get("nt1"))
        p2 = dssr_id_to_text(pair.get("nt2"))
        pairs_dict[p1] = p2
        pairs_dict[p2] = p1
    
    # add base pairing, self edges, and backbone edges
    for i,nt in enumerate(dssr["nts"]):
        nt_id = dssr_id_to_text(nt.get("nt_id"))
        if nt_id in pairs_dict:
            # base pairing edge
            pair_edge = (nt_id, pairs_dict[nt_id])
        else:
            # self edge
            pair_edge = (nt_id, nt_id)
        pairs.append(pair_edge)

        #backbone edges
        if i != (num_nts-1): # if not last nucleotide
            nt_next = dssr_id_to_text(dssr["nts"][i+1].get("nt_id"))
            c1 = nt_id.split(":")[0]
            c2 = nt_next.split(":")[0]
            if c1 == c2: # on same chain
                backbone_edges.append((nt_id, nt_next))
        
    interaction_edges = []
    for key, val in protein_interactions.items():
        for v in val:
            nt = ":".join((key.split(":")[:-1])) #be careful here
            # ss = v.split(":")[3]
            interaction_edges.append((nt, v))
            
    interaction_edges = list(set(interaction_edges))
    # print(interaction_edges) #('C:C', 'A:PRO:277:H')
    stacks = get_stacking_interactions(dssr, ss_dict)

    return pairs,backbone_edges, interaction_edges,stacks