from Bio.PDB import NeighborSearch as Nsearch
from pymol import cmd
from pymol import stored
import sys

"""
NOTE: ONLY WORKS ON PDB FILES FOR NOW!
"""
def getProteinSecStructure(protein, prefix):
    # not needed
    #protein.save("./dssr_output/{}.protein.tmp.cif".format(prefix))
    cmd.load('./dssr_output/{}.tmp.cif'.format(prefix), prefix)
    stored.ss = []
    cmd.iterate( '(n. CA)', "stored.ss.append((chain+':'+resi + ':'+ss))")
    

    residue_seq_to_ss = {} # residue seq (1...) : secondary structure
    
    # go through the large string of stored.ss and extract each secondary structure
    for c in stored.ss:
        residue_seq_to_ss[":".join(c.split(":")[:2]) ] = c.split(":")[2]
    return residue_seq_to_ss

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
    try:
        if atomname in MAJOR_GROOVE_ATOMS[resname]:
            return "major"
        elif atomname in MINOR_GROOVE_ATOMS[resname]:
            return "minor"
        elif "'" in atomname or "P" in atomname:
            return "backbone"
        else:
            return "other"
    except:
        return "other"

"""
Returns dict of interactions of nucleotide with pchnaem, presname, presid[1]
"""
def getInteractions(protein, rna, prefix):
    interactions = {}
    patoms = list(protein.get_atoms())
    if len(patoms) == 0:
        return interactions, {} 
    ns = Nsearch(patoms)
    interactions = {} 
    cut_off = 4

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
        
        
        if(resname in ["HOH","CA"]): #TODO hack for lab meeting
            continue
        
        neighbors = ns.search(atom.coord, radius=cut_off, level="R")
        result = []
        int_type = check_interaction_type(resname, atomname)
        for item in neighbors:
            pchname = item.get_parent().id
            presid = item.get_id()
            presname = item.get_resname()
            try:
                ss = secondary_structure_dict[pchname + ":" + str(presid[1])]
            except:
                ss = 'X' #pymol couldn't return
            result.append("{}:{}:{}:{}".format(pchname, presname, presid[1], ss)) #'A:LEU:269:H'
        
        interactions["{}:{}:{}:{}".format(chname, resname, resid[1], int_type)] = result
    return interactions, secondary_structure_dict
