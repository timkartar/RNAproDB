from Bio.PDB import NeighborSearch as Nsearch
#from pymol import cmd
#from pymol import stored
import os, sys, copy
from utilities import nt_colors, chem_components
from Bio.PDB import PDBIO
from Bio.PDB.DSSP import dssp_dict_from_pdb_file

io = PDBIO()
backend =  os.path.dirname(os.path.abspath(__file__))
frontend = backend + "../rnaprodb_frontend/"
"""
NOTE: ONLY WORKS ON PDB FILES FOR NOW!
"""
def getProteinSecStructure(protein, prefix):
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
    chain_map = {}
    protein = copy.deepcopy(protein)
    for chain in protein.get_chains():
        if len(chain.id) == 1:
            chain_map[chain.id] = chain.id
            alphabet = alphabet.replace(chain.id,"")
        else:
            chain_map[alphabet[0]] = chain.id
            chain.id = alphabet[0]
            alphabet = alphabet.replace(chain.id,"")
    
    io.set_structure(protein)
    io.save('{}-protein.pdb'.format(prefix))
    
    residue_seq_to_ss = {}
    try:
        dssp_tuple = dssp_dict_from_pdb_file('{}-protein.pdb'.format(prefix),DSSP="./external/dssp")
        dssp_dict = dssp_tuple[0]
        for key in dssp_dict.keys():
            chid = key[0]
            resnum = str(key[1][1])
            icode = key[1][2]
            residue_seq_to_ss[":".join([chid, resnum])] = dssp_dict[key][1]
        return residue_seq_to_ss
        #exit()
    except:
        return {}
        
    '''
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
    '''

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
    #for res in protein.get_residues():
    #    print(res.get_full_id())
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
        
        
        if(resname in ["HOH","CA","NA","MG","ZN","ATP"]): #TODO hack for lab meeting
            continue
        if resname not in nt_colors.keys() and resname not in chem_components.keys():
            continue
        neighbors = ns.search(atom.coord, radius=cut_off, level="R")
        result = []
        int_type = check_interaction_type(resname, atomname)
        for item in neighbors:
            pchname = item.get_parent().id
            presid = item.get_id()
            presname = item.get_resname()
            #print(secondary_structure_dict)
            #exit()
            try:
                ss = secondary_structure_dict[pchname + ":" + str(presid[1])]
            except:
                ss = 'X' #pymol couldn't return
            result.append("{}:{}:{}:{}".format(pchname, presname, presid[1], ss)) #'A:LEU:269:H'
        
        interactions["{}:{}:{}:{}:{}".format(chname, resname, resid[1], resid[2].replace(" ",""), int_type)] = result
    return interactions, secondary_structure_dict
