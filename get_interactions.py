from Bio.PDB import NeighborSearch as Nsearch

def getInteractions(protein, rna):
    interactions = {}
    ns = Nsearch(list(protein.get_atoms()))

    interactions = {} 
    cut_off = 5
    for atom in rna.get_atoms():
        chname = atom.get_parent().get_parent().id
        resid = atom.get_parent().get_id()
        resname = atom.get_parent().get_resname()
        atomname = atom.name
        
        
        neighbors = ns.search(atom.coord, radius=cut_off, level="R")
        result = []
        for item in neighbors:
            pchname = item.get_parent().id
            presid = item.get_id()
            presname = item.get_resname()
            #result.append("{}:{}{}".format(pchname, presname, presid))
            result.append("{}:{}:{}".format(pchname, presname, presid[1]))
        #interactions["{}:{}{}:{}".format(chname, resname, resid, atomname)] = result
        interactions["{}:{}:{}:{}".format(chname, resname, resid[1], atomname)] = result
    return interactions


