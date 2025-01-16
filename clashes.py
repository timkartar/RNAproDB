#!/usr/bin/env python3

import os, sys
from Bio.PDB import MMCIFParser, NeighborSearch

# if len(sys.argv) < 2:
#     print("Please provide a PDB ID.")
#     exit(0)

# pdb_id = sys.argv[1]

# curr_dir = os.path.dirname(os.path.abspath(__file__))
# cifs_dir = os.path.dirname(curr_dir) + '/cifs/'
# cif_file = cifs_dir + f'{pdb_id.upper()}.cif'

# parser = MMCIFParser(QUIET=True)

# structure = parser.get_structure(pdb_id, os.path.join(cifs_dir, cif_file))

# van der Waals radius
VDW_RADII = {
    "H": 1.10, "C": 1.70, "N": 1.55, "O": 1.52, "P": 1.80, "S": 1.80, "F": 1.70, "CL": 1.75, "MG": 1.73, "ZN": 1.39, 
    "K": 2.75, "SE": 1.90, "CA": 2.31, "I": 1.98, "U": 1.86, "NA": 2.27, 
}

# ref
# https://www.cgl.ucsf.edu/chimerax/docs/user/commands/clashes.html
# https://www.cgl.ucsf.edu/chimera/docs/ContributedSoftware/findclash/findclash.html
def is_clash(atom1, atom2, overlap_cutoff=0.6, hbond_allowance=0.4):
    dist = atom1 - atom2

    # default radii 1.50
    r1 = VDW_RADII.get(atom1.element, 1.50)
    r2 = VDW_RADII.get(atom2.element, 1.50)
    sum_of_radii = r1 + r2

    # overlap = r1 + r2 - dist
    # clashes criterion: overlap_cutoff < overlap => dist < r1 + r2 - cutoff - hbond allowance
    return dist < (sum_of_radii - overlap_cutoff - hbond_allowance)

def is_donAcc(hbonds, atom1, atom2):
    """Check if atom1 and atom2 form a donor-acceptor pair"""
    atom_pair = (atom1.serial_number, atom2.serial_number)
    return atom_pair in hbonds # and hbonds[atom_pair] != 'questionable'

def get_hbonds(data) -> list:
    hbond_data  = {}
    if 'hbonds' not in data:
        return hbond_data

    # group atom pairs by serial number
    for bond in data['hbonds']:
        atom1 = bond["atom1_serNum"]
        atom2 = bond["atom2_serNum"]
        donAcc_type = bond["donAcc_type"]

        pair_info = (atom1, atom2)
        pair_info_rev = (atom2, atom1)

        if pair_info not in hbond_data:
            hbond_data[pair_info] = donAcc_type
            hbond_data[pair_info_rev] = donAcc_type

    return hbond_data


def find_clashes(structure, data):
    all_atoms = list(structure.get_atoms())
    hbonds = get_hbonds(data)
    ns = NeighborSearch(all_atoms)

    potential_clashes = set()
    for atom in all_atoms:
        neighbors = ns.search(atom.coord, 4.0)

        for neighbor in neighbors:
            if neighbor is atom:
                continue

            # allow atom and its neighbor to be closer if they are a donor-acceptor pair
            hbond_allowance = 0.0
            if is_donAcc(hbonds, atom, neighbor):
                hbond_allowance = 0.4

            # check if they are potential clash by distance and VDW radii
            if is_clash(atom, neighbor, hbond_allowance=hbond_allowance):
                # Sort atoms by ID to avoid duplicates
                if atom.get_id() < neighbor.get_id():
                    potential_clashes.add((atom, neighbor))
                else:
                    potential_clashes.add((neighbor, atom))

    # ref 2
    # https://www.blopig.com/blog/2023/05/checking-your-pdb-file-for-clashing-atoms/
    clashes = []
    potential_clashes = list(potential_clashes)
    for atom_1, atom_2 in potential_clashes:
        # Exclude clashes from atoms in the same residue or 1 residue apart (connected by backbone)
        residue_1, residue_2 = atom_1.parent, atom_2.parent
        chain_1, chain_2 = residue_1.parent, residue_2.parent

        if chain_1 == chain_2 and abs(residue_1.id[1] - residue_2.id[1]) <= 1:
            continue

        if atom_1 - atom_2 == 0:
            continue

        clashes.append((atom_1, atom_2))


    # Print the results
    # for a1, a2 in clashes:
    #     dist = a1 - a2
    #     print(f"Clash: {a1.fullname} (res {a1.parent.id[1]}) - "
    #         f"{a2.fullname} (res {a2.parent.id[1]}), distance: {dist:.2f} Å")

    # Table Data: Node 1, Node 2, distance, Atom 1, Atom 2
    # table = []
    # table.append(f"Node 1,Node 2,Distance,Atom1,Atom2")
    # for a1, a2 in clashes:
    #     dist = a1 - a2
    #     node1 = f"{a1.parent.parent.id}:{a1.parent.id[1]}:"
    #     node2 = f"{a2.parent.parent.id}:{a2.parent.id[1]}:"
    #     dist_str = f"{dist:.3f} Å"
    #     atom1 = f"{a1.fullname}@{a1.parent.parent.parent.serial_num}..{a1.parent.parent.id}.{a1.parent.resname}.{a1.parent.id[1]}."
    #     atom2 = f"{a2.fullname}@{a2.parent.parent.parent.serial_num}..{a2.parent.parent.id}.{a2.parent.resname}.{a2.parent.id[1]}."
        
    #     # entry = f"{node1}\t{node2}\t{dist_str}  \t{atom1} \t{atom2}"
    #     entry = ','.join([node1, node2, dist_str, atom1, atom2])
    #     table.append(entry)

    return clashes

