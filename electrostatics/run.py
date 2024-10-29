from Bio import PDB
import numpy as np
import os, sys

pdb = sys.argv[1]
dirname = home =  os.path.dirname(os.path.abspath(__file__))

# Initialize the parser
parser = PDB.MMCIFParser(QUIET=True)

# Parse the structure file
structure = parser.get_structure('protein', '/home/db/rnaprodb/output/cifs/{}-assembly1.cif'.format(pdb))
alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
chain_map = {}
for chain in structure.get_chains():
    if len(chain.id) == 1:
        chain_map[chain.id] = chain.id
        alphabet = alphabet.replace(chain.id,"")
    else:
        chain_map[alphabet[0]] = chain.id
        chain.id = alphabet[0]
        alphabet = alphabet.replace(chain.id,"")

# Initialize the output PDBIO class
io = PDB.PDBIO()

# Define a custom Select subclass to filter out solvent and ions
class NonSolventSelect(PDB.Select):
    def accept_residue(self, residue):
        # Exclude water molecules and common ions
        solvent_names = ['HOH']  # For water
        ion_names = ['NA', 'K', 'CL', 'MG', 'CA']  # Add any other ions you want to filter out

        # Return True for non-solvent, non-ion residues
        return residue.get_resname() not in solvent_names + ion_names

# Save only the filtered structure (non-solvent, non-ion residues)
io.set_structure(structure)
io.save(dirname + '/full_{}.pdb'.format(pdb), select=NonSolventSelect())
chem_components = dict(np.load("/srv/www/rnaprodb/rnaprodb_dev" + "/modified_parents.npz",allow_pickle=True))
chem_components["AMP"] = "A"
chem_components["CMP"] = "C"
chem_components["GMP"] = "G"
chem_components["TMP"] = "T"
chem_components["UMP"] = "U"

nts = "A,C,G,U,T,DA,DC,DG,DT,DU".split(",")

parser = PDB.PDBParser()
structure = parser.get_structure('protein', 'full_{}.pdb'.format(pdb))

# Define a custom Select subclass to filter out solvent and ions
class NASelect(PDB.Select):
    def accept_residue(self, residue):
        # Exclude water molecules and common ions
        solvent_names = ['HOH']  # For water
        ion_names = ['NA', 'K', 'CL', 'MG', 'CA']  # Add any other ions you want to filter out
        proteins = []
        # Return True for non-solvent, non-ion residues
        name = residue.get_resname()
        if name in nts:
            return True
        else:
            if name in chem_components and chem_components[name] in nts:
                return True
            else:
                return False

io.set_structure(structure)
io.save(dirname + '/na_{}.pdb'.format(pdb), select=NASelect())

if os.path.getsize(dirname + '/full_{}.pdb'.format(pdb)) == os.path.getsize(dirname + '/na_{}.pdb'.format(pdb)):
    os.remove(dirname + '/na_{}.pdb'.format(pdb))
    sys.exit()

class PSelect(PDB.Select):
    def accept_residue(self, residue):
        # Exclude water molecules and common ions
        solvent_names = ['HOH']  # For water
        ion_names = ['NA', 'K', 'CL', 'MG', 'CA']  # Add any other ions you want to filter out
        proteins = []
        # Return True for non-solvent, non-ion residues
        name = residue.get_resname()
        if name in nts:
            return False
        else:
            if name in chem_components and chem_components[name] in nts:
                return False
            else:
                return True


io.save(dirname + '/pro_{}.pdb'.format(pdb), select=PSelect())


