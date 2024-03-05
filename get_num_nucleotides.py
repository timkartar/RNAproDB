from Bio.PDB.MMCIFParser import MMCIFParser

def test_open(cif_file):
    parser = MMCIFParser()
    structure = parser.get_structure('structure', cif_file)
    return

# 
def count_nucleotides_slow(cif_file):
    parser = MMCIFParser()
    structure = parser.get_structure('structure', cif_file)

    nucleotide_count = 0

    # Iterate through all models, chains, and residues
    for model in structure:
        for chain in model:
            for residue in chain:
                # Check if the residue is a nucleotide
                if residue.get_resname() in ["A", "G", "C", "U", "T"]:
                    nucleotide_count += 1

    return nucleotide_count

def count_nucleotides_fast(cif_file):
    nucleotide_residues = set(["A", "G", "C", "U", "T", "DA", "DG", "DC", "DT"])
    nucleotide_count = 0
    has_metadata = False

    with open(cif_file, 'r') as file:
        for line in file:
            # Look for lines that define residue names
            if line.startswith("_pdbx_poly_seq_scheme.hetero "):
                has_metadata = True
                while True:
                    line = next(file)
                    if line.startswith("#"):  # End of section
                        break
                    residue_name = line.strip().split()[3].strip()
                    if residue_name in nucleotide_residues:
                        nucleotide_count += 1
    if not has_metadata:
        nucleotide_count = count_nucleotides_slow(cif_file)
    return nucleotide_count

# about 30 seconds!

# Replace 'your_file.cif' with the path to your CIF file
file_path = '/home/aricohen/Desktop/django-react-rnaprodb/frontend/public/cifs/6gv4-assembly1.cif'
nucleotides = count_nucleotides_fast(file_path)