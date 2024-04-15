import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from ase.io import read

def plot_energy_from_xyz(xyz_file):
    # Read the XYZ file
    atoms = read(xyz_file, index=':')
    
    # Extract energies and geometry indices
    energies = []
    geometry_indices = []
    
    for i, atom in enumerate(atoms):
        energy = atom.info['energy']
        energies.append(energy)
        geometry_indices.append(i)
    
    # Convert lists to numpy arrays for plotting
    energies = np.array(energies)
    geometry_indices = np.array(geometry_indices)
    
    # Plot energies
    plt.plot(
             geometry_indices,
             energies,
             marker='o',
             markersize=3,
             linewidth=0.5,
             linestyle='-',
             label=os.path.basename(xyz_file))

# Check if two arguments are provided
if len(sys.argv) != 3:
    print("Usage: python script.py <xyz_file1> <xyz_file2>")
    sys.exit(1)

# Plot energies from the first XYZ file
plot_energy_from_xyz(sys.argv[1])

# Plot energies from the second XYZ file
plot_energy_from_xyz(sys.argv[2])

# Set plot labels and title
plt.xlabel('Geometry Index')
plt.ylabel('Energy (eV)')
plt.title('Energy vs. Geometry Index')
plt.legend()
plt.grid(True)

# Save the plot in the same directory as the first XYZ file
output_file = os.path.splitext(sys.argv[1])[0] + '_energy_plot.png'
plt.savefig(output_file)

# Show the plot
plt.show()

