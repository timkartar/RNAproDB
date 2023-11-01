import os, sys, json
import subprocess

def getSS(prefix):
    subprocess.run(["x3dna-dssr","-i=../frontend/public/cifs/{}-assembly1.cif".format(prefix),"--nested"])
    l = open("dssr-2ndstrs.dbn","r").readlines()
    seq = l[1].strip().split("&")
    sec = l[2].strip().split("&")
    return seq, sec

# Return JSON of ss and some metadata
def processSS(ss):
    # first tuple is sequence/s
    num_chains = 0
    ss_json_list = []
    ss_json_object = {}

    for i, chain in enumerate(ss[0]):
        # ignore chain if one resiude, subject to change!
        if len(chain) < 2:
            continue
        cur_ss = {}
        cur_ss['sequence'] = chain
        cur_ss['structure'] = ss[1][i] 
        num_chains += 1
        ss_json_list.append(cur_ss)
    # second tuple is ss
    ss_json_object["ssList"] = ss_json_list
    ss_json_object['numChains'] = num_chains
    # print(ss_json_object)
    return ss_json_object