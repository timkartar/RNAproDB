import re

def extract_first_letter_and_last_number(s):
    match = re.search(r'([A-Z])\..*\.([0-9]+)\.', s)
    if match:
        return f"{match.group(1)}:{match.group(2)}"
    return None


def getLW(dssr):
    lw_dict = {}
    for pair in dssr['pairs']:
        nt1_result = extract_first_letter_and_last_number(pair['nt1'])
        nt2_result = extract_first_letter_and_last_number(pair['nt2'])
        lw_dict[(nt1_result, nt2_result)] = pair['LW']
    return lw_dict





