from data import data as D
from Bio.PDB import is_aa
import numpy as np
import os
import copy
from utilities import d3to1, chem_components


def splitEntities(structure, regexes=None, atom_mapper=None, mi=0):
    """Docstring"""
    if regexes is None:
        regexes = D.regexes
    RNA = ["A","C","G","U","T","DA","DC","DG","DT"]    
    pro = []
    lig = []
    for chain in structure[mi]:
        for residue in chain:
            resname = residue.get_resname().strip()
            if resname in d3to1.keys():
                pro.append(residue.get_full_id())
            else:
                if resname in RNA:
                    lig.append(residue.get_full_id())
                elif resname in chem_components.keys():
                    lig.append(residue.get_full_id())
    
    protein = copy.deepcopy(structure)
    rna = copy.deepcopy(structure)


    for res in structure.get_residues():
        if res.get_full_id() not in pro:
            protein[mi][res.get_parent().id].detach_child(res.get_id())

    for res in structure.get_residues():
        if res.get_full_id() not in lig:
            rna[mi][res.get_parent().id].detach_child(res.get_id())

   
    return protein[0], rna[0]
    #return structure.slice(structure, pro, name='{}_protein'.format(structure.name)), structure.slice(structure, lig, name='{}_ligand'.format(structure.name))
