import os
from calcs import calc
from ase.io import read
from ase.neb import NEBTools
import matplotlib.pyplot as plt
from initialize_path import setup_images

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

nebtools_i = NEBTools(initial_images)
Ef_i, dE_i = nebtools_i.get_barrier(fit=False)

nebtools_f = NEBTools(optimized_images)
Ef_f, dE_f = nebtools_f.get_barrier(fit=False)

output_filename = f'neb_band_{method}_NEBOptimizer_None.png'
output_path = os.path.join(logdir, output_filename)

fig, ax = plt.subplots()
nebtools_i.plot_band(ax=ax)

for image in initial_images:
    ax.plot(image.get_potential_energy(), label='Initial Image', linestyle='--')
for image in optimized_images:
    ax.plot(image.get_potential_energy(), label='Optimized Image', linestyle='-')
ax.set_xlabel('Image Index')
ax.set_ylabel('Potential Energy (eV)')
ax.legend()

plt.savefig(output_path)
plt.close(fig)

