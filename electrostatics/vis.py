import trimesh, sys
import numpy as np
import io, json

data = dict(np.load("./{}.npz".format(sys.argv[1]), allow_pickle=True))
data['points'] = data['points'] - np.mean(data['points'], axis=0)
#for key in data:
#    data[key] = data[key].tolist()

#with open("data.json", "w") as outfile: 
#    json.dump(data, outfile)

mesh = trimesh.Trimesh(vertices=data['points'],
                       faces=data['faces'],process=False)

import matplotlib.pyplot as plt

import matplotlib.cm as cm
cmap = cm.bwr_r
values = data["pdata"]
norm = plt.Normalize(-20, +20)
colors = cmap(norm(values))
mesh.visual.vertex_colors = colors

mesh.export('/home/db/rnaprodb/output/electrostatics/{}.ply'.format(sys.argv[1]))

