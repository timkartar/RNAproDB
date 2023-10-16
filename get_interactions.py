from Bio.PDB import NeighborSearch as Nsearch
from pymol import cmd
from pymol import stored

"""
NOTE: ONLY WORKS ON PDB FILES FOR NOW!
"""
def getProteinSecStructure(protein, prefix):
    # not needed
    #protein.save("./dssr_output/{}.protein.tmp.cif".format(prefix))
    cmd.load('./dssr_output/{}.tmp.cif'.format(prefix), prefix)
    stored.ss = ""
    cmd.iterate( '(n. CA)', 'stored.ss = stored.ss + ("%1s"%ss)')
    counter = 1

    residue_seq_to_ss = {} # residue seq (1...) : secondary structure

    # go through the large string of stored.ss and extract each secondary structure
    for c in stored.ss:
        residue_seq_to_ss[counter] = c
        counter += 1
    # vals = set()
    # for i in range(1,len(residue_seq_to_ss)+1):
    #     if residue_seq_to_ss[i] not in vals:
    #         vals.add(residue_seq_to_ss[i])
    #         print(residue_seq_to_ss[i])
    # print(vals) # only see S, H, and L
    
    #print(residue_seq_to_ss)
    #import sys
    #sys.exit()
    
    return residue_seq_to_ss
    # resseq = residue.get_full_id()[3][1]

def check_interaction_type(resname, atomname):
    '''
    check whether the protein RNA interaction is in major groove/minor groove/backbone
    '''
    MAJOR_GROOVE_ATOMS = {
    "A": ["N7", "N6"],
    "C": ["N4", "C5"],
    "G": ["N7", "O6"],
    "U": ["O4", "C5"]
    }

    MINOR_GROOVE_ATOMS = {
    "A": ["N3", "C2"],
    "C": ["O2"],
    "G": ["N3", "C2"],
    "U": ["O2"]
    }
    if atomname in MAJOR_GROOVE_ATOMS[resname]:
        return "major"
    elif atomname in MINOR_GROOVE_ATOMS[resname]:
        return "minor"
    elif "'" in atomname or "P" in atomname:
        return "backbone"
    else:
        return "other"

"""
Returns dict of interactions of nucleotide with pchnaem, presname, presid[1]
"""
def getInteractions(protein, rna, prefix):
    interactions = {}
    ns = Nsearch(list(protein.get_atoms()))
    interactions = {} 
    cut_off = 5

    secondary_structure_dict = getProteinSecStructure(protein, prefix) # residue number to secondary structure abbreviation 
    # print("secondary structure length from pymol: {}".format(len(secondary_structure_dict)))
    # num_resiudes_biopython = 0
    # for residue in protein.get_residues():
    #     num_resiudes_biopython += 1
    # print("num resiudes biopython: {}".format(num_resiudes_biopython)) # Lengths match up!
    for atom in rna.get_atoms():
        chname = atom.get_parent().get_parent().id
        resid = atom.get_parent().get_id()
        resname = atom.get_parent().get_resname()
        atomname = atom.name
        
        
        neighbors = ns.search(atom.coord, radius=cut_off, level="R")
        result = []
        int_type = check_interaction_type(resname, atomname)
        for item in neighbors:
            pchname = item.get_parent().id
            presid = item.get_id()
            presname = item.get_resname()
            ss = secondary_structure_dict[presid[1]]
            result.append("{}:{}:{}:{}".format(pchname, presname, presid[1], ss)) #'A:LEU:269:H'
        
        interactions["{}:{}:{}:{}".format(chname, resname, resid[1], int_type)] = result
    return interactions, secondary_structure_dict
