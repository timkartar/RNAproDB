from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.Polypeptide import is_aa
from tqdm import tqdm
import pickle

def aa_counter(cif_file):
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("structure", cif_file)    
    aa_count = 0
    model = structure[0]
    for chain in model:
        for residue in chain:
            if is_aa(residue):
                aa_count += 1
    return aa_count

if __name__ == "__main__":
    with  open('hirad/nakb_prna_ids.txt', 'r') as file:
        PDB_IDS = file.read().split(',')
    aa_counts_dict = {}
    for ID in tqdm(PDB_IDS, desc='Parsing through PDB IDs'):
        try:
            aa_count = aa_counter(f'hirad/CIF_structs/{ID}-assembly1.cif')
            aa_counts_dict[ID] = aa_count
        except Exception as e:
            print(f"An error occurred for PDB {ID}: {e}")
    with open('aa_counts_dict.pickle', 'wb') as file:
        pickle.dump(aa_counts_dict, file)
