import os, sys
import json

home =  os.path.dirname(os.path.abspath(__file__)) #change this line only 
pdb_path = "{}/output/".format(home)


def readJSON(pdbid):
    # with open('{}/output/{}_graph.json'.format(home, pdbid), 'r') as infile:
    #     data_string = infile.read()
    with open("{}/output/{}_graph.json".format(home,pdbid), 'r') as infile:
        json_output = json.load(infile)
        return json_output

def convert_coordinates(data, algorithm):
    if algorithm == "PCA" or algorithm == "None":
        for node in data['nodes']:
            node['x'] = node['pca_x']
            node['y'] = node['pca_y']
    elif algorithm == "RNAScape":
        for node in data['nodes']:
            node['x'] = node['rnascape_x']
            node['y'] = node['rnascape_y']
    elif algorithm == "SecondaryStructure":
        for node in data['nodes']:
            node['x'] = node['viennarna_x']
            node['y'] = node['viennarna_y']
    return data


if __name__ == "__main__":
    data = readJSON(sys.argv[1]) # assumes pre-computed JSON file, reads in pdbid
    algorithm = sys.argv[2]
    data = convert_coordinates(data, algorithm)
    final_json_str = json.dumps(data)
    print(final_json_str) #use console output
