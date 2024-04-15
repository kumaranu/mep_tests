import os
from calcs import calc
from ase.io import read
from ase.neb import ElasticBand
import matplotlib.pyplot as plt
from initialize_path import setup_images
import numpy as np

method = 'aseneb'
logdir = '/global/home/users/kumaranu/Documents/gpu_jobs'

r_p = read(os.path.join(logdir, 'r_p.xyz'), ':')
intermediate_images = read(os.path.join(logdir, 'output.xyz'), ':')
optimized_images = read(os.path.join(logdir, 'optimized_path.xyz'), ':')

reactant = r_p[0]
reactant.calc = calc()
product = r_p[1]
product.calc = calc()

initial_images = [reactant]
for image in intermediate_images:
    image.calc = calc()
    initial_images.append(image)
initial_images.append(product)

# Calculate RMSD values for initial images
initial_rmsd = []
for image in initial_images:
    initial_rmsd.append(image.get_all_distances(product).sum() / len(reactant))

# Calculate RMSD values for optimized images
optimized_rmsd = []
for image in optimized_images:
    optimized_rmsd.append(image.get_all_distances(product).sum() / len(reactant))

output_filename = f'neb_band_{method}_NEBOptimizer_None.png'
output_path = os.path.join(logdir, output_filename)

fig, ax = plt.subplots()

# Plot initial images
ax.plot(initial_rmsd, [image.get_potential_energy() for image in initial_images], marker='o', linestyle='-', label='Initial Image', color='blue')

# Plot optimized images on top of initial images
ax.plot(optimized_rmsd, [image.get_potential_energy() for image in optimized_images], marker='o', linestyle='-', label='Optimized Image', color='red')

ax.set_xlabel('RMSD from Reactant')
ax.set_ylabel('Potential Energy (eV)')
ax.legend()

plt.savefig(output_path)
plt.close(fig)

