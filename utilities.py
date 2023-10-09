"""
Input: ('p'/'nt', name, position, chain, ss (protein only))
Returns text string from node tuple representation
"""
def node_to_text(node):
    if len(node) == 5: # IS A PROTEIN
        return node[3] + ":" + node[1] + ":" + node[2] + ":" + node[4]
    return node[3] + ":" + node[1] + ":" + node[2] # is a NT