from data import data
#### Biopython Disordered Atom Fix ####
try:
    import Bio.PDB
    copy = Bio.PDB.Atom.copy
    def myCopy(self):
        shallow = copy.copy(self)
        for child in self.child_dict.values():
            shallow.disordered_add(child.copy())
        return shallow
    Bio.PDB.Atom.DisorderedAtom.copy=myCopy
except ModuleNotFoundError:
    # BioPython isn't found - futher attempts to import will raise an exception
    pass


def isNucleotide(resname):
    if resname in data.standard_RNA_nucleotides:
        # it's standard, not modified
        return True
    
    if resname in data.chem_components and '_chem_comp.mon_nstd_parent_comp_id' in data.chem_components[resname]:
        return (data.chem_components[resname]['_chem_comp.mon_nstd_parent_comp_id'] in data.standard_RNA_nucleotides)
    else:
        # has no standard parent field - can't be modified
        return False

def cleanRNA(structure,
        fix_modified_nucleotide_hetflags=False,
        remove_hetatoms=True
    ):
    for chain in structure.get_chains():
        modified = []
        pos = 0
        modified_positions = []
        remove = []
        for residue in chain:
            resname = residue.get_resname()
            if isNucleotide(resname) and residue.get_id()[0][0] == "H":
                modified.append(residue)
                modified_positions.append(pos)
            pos += 1
            if not isNucleotide(resname):
                remove.append(residue)
        if fix_modified_nucleotide_hetflags:
            # remove heteratom flags on modified nucleotides
            pos_idx = 0
            for residue in modified:
                rid = residue.get_id()
                chain.detach_child(rid)
                residue.id = (' ', rid[1], rid[2])
                chain.insert(modified_positions[pos_idx], residue)
                
                if residue.get_resname() not in data.standard_RNA_nucleotides:
                    residue.resname = data.chem_components[residue.get_resname()]['_chem_comp.mon_nstd_parent_comp_id']
                
                pos_idx += 1
        
        
        for residue in remove:
            rid = residue.get_id()
            chain.detach_child(rid)

        if remove_hetatoms:
            for residue in modified:
                rid = residue.get_id()
                chain.detach_child(rid)
    
    toremove_chains = []
    for chain in structure.get_chains():
        if len(chain.child_list) == 0:
            toremove_chains.append(chain.id)

    for item in toremove_chains:
        structure.structure[0].detach_child(item)
    return structure
