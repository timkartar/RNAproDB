import re
from utilities import dssr_id_to_text
def extract_first_letter_and_last_number(s):
    #match = re.search(r'([A-Z])\..*\.([0-9]+)\.', s)
    #if match:
    #    return f"{match.group(1)}:{match.group(2)}"
    #return None
    spl = dssr_id_to_text(s).split(":")
    return ":".join([spl[0],spl[2],spl[3]])


def getLW(dssr):
    lw_dict = {}
    if 'pairs' not in dssr.keys():
        return lw_dict
    for pair in dssr['pairs']:
        nt1_result = extract_first_letter_and_last_number(pair['nt1'])
        nt2_result = extract_first_letter_and_last_number(pair['nt2'])
        lw_dict[(nt1_result, nt2_result)] = pair['LW']
    #print(lw_dict)
    return lw_dict





