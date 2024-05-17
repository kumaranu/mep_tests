import glob, os, sys
from ase.io import read
from calcs import calc
import numpy as np
from ase import Atoms
from multiprocessing import Pool


def align_geom(refgeom, geom):
    """Find translation/rotation that moves a given geometry to maximally overlap
    with a reference geometry. Implemented with Kabsch algorithm.

    Args:
        refgeom:    The reference geometry to be rotated to
        geom:       The geometry to be rotated and shifted

    Returns:
        RMSD:       Root-mean-squared difference between the rotated geometry
                    and the reference
        new_geom:   The rotated geometry that maximumally overal with the reference
    """
    center = np.mean(refgeom, axis=0)   # Find the geometric center
    ref2 = refgeom - center
    geom2 = geom - np.mean(geom, axis=0)
    cov = np.dot(geom2.T, ref2)
    v, sv, w = np.linalg.svd(cov)

    if np.linalg.det(v) * np.linalg.det(w) < 0:
        sv[-1] = -sv[-1]
        v[:, -1] = -v[:, -1]
    u = np.dot(v, w)
    new_geom = np.dot(geom2, u) + center
    rmsd = np.sqrt(np.mean((new_geom - refgeom) ** 2))
    return rmsd, new_geom


def process_file(i):
    try:
        if os.path.exists(f'{i:03}'):
            os.chdir(f'{i:03}')
            if os.path.exists('neb_band_aseneb_NEBOptimizer_None.txt'):
                ase_traj = read('optimized_path_aseneb_NEBOptimizer_None.xyz', ':')
                energies = [i.get_potential_energy() for i in ase_traj]
                ts_index = np.argmax(energies)
                
                mlcalculator = calc()
                ts_atoms = read('TS.xyz')
                mlcalculator.calculate(ts_atoms)
                energy_ts = mlcalculator.results['energy']
                rmsd, new_geom = align_geom(ase_traj[ts_index].get_positions(), ts_atoms.get_positions())
                
                fwd_barrier_neb = energies[ts_index] - energies[0]
                rev_barrier_neb = energies[ts_index] - energies[-1]
                ref_fwd_barrier = energy_ts - energies[0]
                ref_rev_barrier = energy_ts - energies[-1]

                ts_e_err = energy_ts - energies[ts_index]
                fwd_barrier_err = fwd_barrier_neb - ref_fwd_barrier
                rev_barrier_err = rev_barrier_neb - ref_rev_barrier
                os.chdir('../')
                return f'{i:03}: TS energy err: {ts_e_err:.6f} fwd barrier err: {fwd_barrier_err:.6f} Rev barrier err: {rev_barrier_err:.6f} Geometry RMSD: {rmsd:.6f}'
            else:
                os.chdir('../')
    except Exception as e:
        return f'Error processing file {i}: {str(e)}'


if __name__ == '__main__':

    logdir = '/global/scratch/users/kumaranu/saved_files/neb_nn_inputs'
    os.chdir(logdir)

    with Pool(5) as pool:
        results = pool.map(process_file, range(265))
   
    missing_indices = [i for i, result in enumerate(results) if result is None]
    if missing_indices:
        print(f'Missing results for indices: {missing_indices}')
 
    for result in results:
        if result is not None:
            print(result)
