import os, sys
lines = open("./components.cif","r").readlines()

res = []
parent = []
d = {}
for l in lines:
    if "_chem_comp.id" in l:
        res.append(l.strip().split(" ")[-1])
    if "_chem_comp.mon_nstd_parent_comp_id" in l:
        parent.append(l.strip().split(" ")[-1])
for item in range(len(res)):
    if parent[item] != "?":
        d[res[item]] = parent[item]

import numpy as np
np.savez("modified_parents.npz",**d)
