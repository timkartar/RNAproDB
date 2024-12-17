import trimesh
import numpy as np
import matplotlib.pyplot as plt
import sys

# Load the .ply file
input_ply_path = sys.argv[1]
output_ply_path = sys.argv[2]

# Load the mesh
mesh = trimesh.load(input_ply_path)

# Check if the mesh has vertex colors
#if 'vertex_colors' not in mesh.metadata:
#    print("The PLY file does not contain color information.")
#else:
# Extract RGB colors from vertex colors (assuming colors are in the format [R, G, B, A])
colors = mesh.visual.vertex_colors[:, :3] / 255.0  # Normalize to range [0, 1]

# Convert colors to grayscale
grayscale = np.dot(colors, [0.2989, 0.587, 0.114])  # Luminosity method for grayscale

# Apply a new colormap (e.g., 'viridis')
colormap = plt.cm.get_cmap('bwr_r')
#from matplotlib.colors import LinearSegmentedColormap

# Create a custom colormap
#colors = [(0, 'red'), (0.2, 'red'), (0.4, 'white'), (0.6, 'white'), (0.8, 'blue'), (1, 'blue')]
#colormap = LinearSegmentedColormap.from_list('custom_cmap', colors)

new_colors = colormap(grayscale)[:, :3] * 255  # Apply colormap and convert to 0-255 range

# Update the mesh's vertex colors
mesh.visual.vertex_colors[:, :3] = new_colors.astype(np.uint8)

# Export the modified mesh to a new .ply file
mesh.export(output_ply_path)
print(f"New PLY file with remapped colors saved as '{output_ply_path}'.")

