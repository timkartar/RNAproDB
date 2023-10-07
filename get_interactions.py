from Bio.PDB import NeighborSearch as Nsearch
from pymol import cmd
from pymol import stored
from pymol import selector
from Bio.PDB import MMCIFParser
from Bio.PDB import PDBParser

"""
NOTE: ONLY WORKS ON PDB FILES FOR NOW!
"""
def getProteinSecStructure(protein, prefix):
    protein.save_pdb("./dssr_output/{}.tmp.pdb".format(prefix))
    cmd.load('./dssr_output/{}.tmp.pdb'.format(prefix), prefix)
    stored.ss = ""
    cmd.iterate( '(n. CA)', 'stored.ss = stored.ss + ("%1s"%ss)')
    counter = 1

    residue_seq_to_ss = {} # residue seq (1...) : secondary structure

    for c in stored.ss:
        residue_seq_to_ss[counter] = c
        counter += 1
    return residue_seq_to_ss
    # resseq = residue.get_full_id()[3][1]


"""
Returns dict of interactions of nucleotide with pchnaem, presname, presid[1]
"""
def getInteractions(protein, rna, prefix):
    interactions = {}
    ns = Nsearch(list(protein.get_atoms()))
    interactions = {} 
    cut_off = 5

    secondary_structure_dict = getProteinSecStructure(protein, prefix) # residue number to secondary structure abbreviation 


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
            result.append("{}:{}:{}:{}".format(pchname, presname, presid[1], secondary_structure_dict[presid[1]])) #'A:LEU:269:H'
            #return secondary structure
        #interactions["{}:{}{}:{}".format(chname, resname, resid, atomname)] = result
        interactions["{}:{}:{}:{}".format(chname, resname, resid[1], atomname)] = result
    return interactions