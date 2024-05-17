import os
import sys
import numpy as np
from calcs import calc
from ase.io import read
import matplotlib.pyplot as plt


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

    # Find the index and energy value of the maximum energy point
    max_energy_idx = np.argmax(energies)
    max_energy = energies[max_energy_idx]

    # Plot energies
    plt.plot(
             geometry_indices,
             energies,
             marker='o',
             markersize=3,
             linewidth=0.5,
             linestyle='-',
             label=os.path.basename(xyz_file),
    )
    return energies[0], max_energy, energies[-1]


def gen_energy_profile(
        geodesic_xyz,
        neb_path_xyz,
        ts_xyz,
):
    # Plot energies from the first XYZ file
    _, _, _ = plot_energy_from_xyz(geodesic_xyz)

    # Plot energies from the second XYZ file
    e_r, e_highest, e_p = plot_energy_from_xyz(neb_path_xyz)

    ts = read(ts_xyz)
    mlcalculator = calc()
    mlcalculator.calculate(ts)
    e_ts = mlcalculator.results['energy']

    # Set plot labels and title
    plt.xlabel('Geometry Index')
    plt.ylabel('Energy (eV)')
    plt.title(f'neb barrier f: {e_highest-e_r:.3f} eV, '
              f'Ref. barrier f: {e_ts-e_r:.3f} eV\n'
              f'neb barrier r: {e_highest-e_p:.3f} eV, '
              f'Ref. barrier r: {e_ts-e_p:.3f} eV')
    plt.legend(['geodesic', 'neb'])
    plt.grid(True)

    # Save the plot in the same directory as the first XYZ file
    output_file = os.path.splitext(geodesic_xyz)[0] + '_energy_plot.png'
    plt.savefig(output_file)

    # Show the plot
    plt.show()
    return e_highest-e_ts


if __name__ == '__main__':
    # Check if two arguments are provided
    # if len(sys.argv) != 3:
    #     print("Usage: python script.py <xyz_file1> <xyz_file2>")
    #     sys.exit(1)
    #ref_dir = '/home/kumaranu/Downloads/neb_nn_inputs'
    #index = '001'
    ref_dir = '/home/kumaranu/Downloads/neb_rgd1_inputs'
    for i in range(999):
        index = f'{i:03}'
        # index = '679'
        try:
            geodesic_xyz = os.path.join(
                ref_dir,
                index,
                'geodesic_path.xyz',
            )
            neb_path_xyz = os.path.join(
                ref_dir,
                index,
                'optimized_path_aseneb_NEBOptimizer_None.xyz'
            )
            ts_xyz = os.path.join(
                ref_dir,
                index,
                'TS.xyz',
            )

            barrier_err = gen_energy_profile(
                geodesic_xyz,
                neb_path_xyz,
                ts_xyz,
            )
        except:
            continue
