import os, sys, json
import subprocess
from utilities import dssr_id_to_text
def getSS(prefix, data):
    #subprocess.run(["x3dna-dssr","-i=../rnaprodb_frontend/public/cifs/{}-assembly1.cif".format(prefix),"--nested"])
    #try:
    #print(json.dumps(data, indent=4))
    l = ["",data["dbn"]["all_chains"]["bseq"], data["dbn"]["all_chains"]["sstr"]]
    seq = l[1].strip().split("&")
    sec = l[2].strip().split("&")
    #except:
    #    return [],[]
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

def dssrSS(data):
    dssrss = {}
    if 'nts' not in data.keys():
        return dssrss
    else:
        for item in data['nts']:
            id1 = dssr_id_to_text(item['nt_id'])
            dssrss[id1] = []

    try:
        if 'stems' in data.keys():
            for item in data['stems']:
                stem_index = item['index']
                helix_index = item['helix_index']
                for pair in item['pairs']:
                    id1 = dssr_id_to_text(pair['nt1'])
                    id2 = dssr_id_to_text(pair['nt2'])
                    dssrss[id1].append(["stem", stem_index, helix_index])
                    dssrss[id2].append(["stem", stem_index, helix_index])
        if 'hairpins' in data.keys():
            for item in data['hairpins']:
                index = item['index']
                for nt in item['nts_long'].split(","):
                    id1 = dssr_id_to_text(nt)
                    dssrss[id1].append(["hairpins", index, item['type']])
        if 'junctions' in data.keys():
            for item in data['junctions']:
                index = item['index']
                for nt in item['nts_long'].split(","):
                    id1 = dssr_id_to_text(nt)
                    dssrss[id1].append(["junctions", index, item['type']])
        if 'bulges' in data.keys():
            for item in data['bulges']:
                index = item['index']
                for nt in item['nts_long'].split(","):
                    id1 = dssr_id_to_text(nt)
                    dssrss[id1].append(["bulges", index, item['type']])
        
        if 'iloops' in data.keys():
            for item in data['iloops']:
                index = item['index']
                for nt in item['nts_long'].split(","):
                    id1 = dssr_id_to_text(nt)
                    dssrss[id1].append(["iloops", index, item['type']])
        if 'ssSegments' in data.keys():
            for item in data['ssSegments']:
                index = item['index']
                for nt in item['nts_long'].split(","):
                    id1 = dssr_id_to_text(nt)
                    dssrss[id1].append(["ssSegments", index, 'ssSegments'])
    except:
        pass
    ''' Tertiary structure
    if 'Kturns' in data.keys():
        for item in data['kturns']:
            index = item['index']
    '''
    return dssrss
