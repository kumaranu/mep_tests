import os
from calcs import calc
from ase.io import read, write
from sella_wrapper import sella_wrapper


def setup_images(
        logdir: str,
        xyz_r_p: str,
        n_intermediate: int = 40,
):
    """
    Sets up intermediate images for NEB calculations between reactant and product states.

    Parameters:
    logdir (str): Directory to save the intermediate files.
    xyz_r_p (str): Path to the XYZ file containing reactant and product structures.
    n_intermediate (int): Number of intermediate images to generate.

    Returns:
    List: List of ASE Atoms objects with calculated energies and forces.
    """
    try:
        # Ensure the log directory exists
        os.makedirs(logdir, exist_ok=True)

        # Read reactant and product structures
        reactant = read(xyz_r_p, index='0')
        product = read(xyz_r_p, index='1')

        # Optimize reactant and product structures using sella
        for atom, name in zip([reactant, product], ['reactant', 'product']):
            atom.calc = calc()
            traj_file = os.path.join(logdir, f'{name}_opt.traj')
            sella_wrapper(atom, traj_file=traj_file, sella_order=0)

        # Save optimized reactant and product structures
        r_p_path = os.path.join(logdir, "r_p.xyz")
        write(r_p_path, [reactant.copy(), product.copy()])

        # Generate intermediate images using geodesic interpolation
        output_xyz = os.path.join(logdir, 'output.xyz')
        os.system(f'geodesic_interpolate {r_p_path} --output {output_xyz} --nimages {n_intermediate}')
        images = read(output_xyz, index=':')

        # Calculate energies and forces for each intermediate image
        for image in images:
            image.calc = calc()
            ml_calculator = calc()
            ml_calculator.calculate(image)
 
            energy = ml_calculator.results['energy']
            forces = ml_calculator.results['forces']
 
            image.info['energy'] = energy
            image.arrays['forces'] = forces
 
        # Save the geodesic path
        geodesic_path = os.path.join(logdir, 'geodesic_path.xyz')
        write(geodesic_path, images)

        return images

    except Exception as e:
        print(f"An error occurred: {e}")
        return []


if __name__ == '__main__':

    setup_images()
