import os
from calcs import calc
from ase import Atoms
from ase.io import read, write
from sella_wrapper import sella_wrapper

from geodesic_interpolate.fileio import read_xyz, write_xyz
from geodesic_interpolate.interpolation import redistribute
from geodesic_interpolate.geodesic import Geodesic


def geodesic_interpolate_wrapper(
    r_p_atoms: Atoms,
    nimages: int = 17,
    sweep: bool = None,
    output: str = "interpolated.xyz",
    tol: float = 2e-3,
    maxiter: int = 15,
    microiter: int = 20,
    scaling: float = 1.7,
    friction: float = 1e-2,
    dist_cutoff: float = 3,
    save_raw: str = None):
    """
    Interpolates between two geometries and optimizes the path.

    Parameters:
    filename (str): XYZ file containing geometries.
    nimages (int): Number of images. Default is 17.
    sweep (bool): Sweep across the path optimizing one image at a time.
                  Default is to perform sweeping updates if there are more than 35 atoms.
    output (str): Output filename. Default is "interpolated.xyz".
    tol (float): Convergence tolerance. Default is 2e-3.
    maxiter (int): Maximum number of minimization iterations. Default is 15.
    microiter (int): Maximum number of micro iterations for sweeping algorithm. Default is 20.
    scaling (float): Exponential parameter for Morse potential. Default is 1.7.
    friction (float): Size of friction term used to prevent very large change of geometry. Default is 1e-2.
    dist_cutoff (float): Cut-off value for the distance between a pair of atoms to be included in the coordinate system. Default is 3.
    save_raw (str): When specified, save the raw path after bisections but before smoothing. Default is None.
    """
    # Read the initial geometries.
    symbols = r_p_atoms[0].get_chemical_symbols()

    X = [conf.get_positions() for conf in r_p_atoms]

    if len(X) < 2:
        raise ValueError("Need at least two initial geometries.")

    # First redistribute number of images. Perform interpolation if too few and subsampling if too many images are given
    raw = redistribute(symbols, X, nimages, tol=tol * 5)
    if save_raw is not None:
        write_xyz(save_raw, symbols, raw)

    # Perform smoothing by minimizing distance in Cartesian coordinates with redundant internal metric
    # to find the appropriate geodesic curve on the hyperspace.
    smoother = Geodesic(symbols, raw, scaling, threshold=dist_cutoff, friction=friction)
    if sweep is None:
        sweep = len(symbols) > 35
    try:
        if sweep:
            smoother.sweep(tol=tol, max_iter=maxiter, micro_iter=microiter)
        else:
            smoother.smooth(tol=tol, max_iter=maxiter)
    finally:
        # Save the smoothed path to output file. try block is to ensure output is saved if one ^C the process, or there is an error
        write_xyz(output, symbols, smoother.path)
    return symbols, smoother.path


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
        symbols, smoother_path =\
            geodesic_interpolate_wrapper([reactant.copy(), product.copy()])
        images = [Atoms(symbols=symbols, positions=conf) for conf in smoother_path]

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
