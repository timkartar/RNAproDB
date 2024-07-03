import json
from make_table import makeTable

with open("/srv/www/rnaprodb/rnaprodb_dev/output/1ivs_pca_graph.json", 'r') as json_file:
        json_output = json.load(json_file)

makeTable(json_output)
