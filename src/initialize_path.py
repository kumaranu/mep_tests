import os
import copy
import numpy as np
from calcs import calc
from ase.io import read, write
from ase.io import Trajectory
from ase.optimize import BFGS, ODE12r
from sella import Sella
from typing import Any, Dict, Iterator, List, Sequence, Tuple, TypeVar, Union
import ase.units as units
from math import pi, sin, sqrt
from ase.units import kcal, mol, Ang
from sella_wrapper import sella_wrapper


def get_hessian(
                ts,
                hess_precision=1e-5,
                method='fwd_diff'
               ):
    n_atoms = len(ts)
    hessian = np.zeros((n_atoms, 3, n_atoms, 3))

    # Calculate Hessian using finite difference
    hess_precision = 1e-5 # Hessian precision
    for A_ in range(n_atoms):
        for X_ in range(3):
            ts.positions[A_][X_] += hess_precision
            ts.set_calculator(calc())
            forces_pos = ts.get_forces()
            ts.positions[A_][X_] -= 2 * hess_precision
            ts.set_calculator(calc())
            forces_neg = ts.get_forces()
            hessian[:, :, A_, X_] = -(forces_pos - forces_neg) / (2 * hess_precision) * (kcal/mol/Ang)
    return hessian


def get_hessian_2d(ts):
    hessian = get_hessian(ts)
    n_atoms = np.shape(hessian)[0]
    return np.reshape(hessian, (3*n_atoms, 3*n_atoms))


def get_eigs(hessian):
    n_atoms = np.shape(hessian)[0]
    eigvals, _ = np.linalg.eig(np.reshape(hessian, (3*n_atoms, 3*n_atoms)))
    eigvals_real = np.real(eigvals)
    formatted_eigvals = ["{:.3e}".format(val) for val in eigvals_real]

    print('eigvals')
    for val in formatted_eigvals:
        print(val)


def _energies_and_modes(ts) -> Tuple[np.ndarray, np.ndarray]:
    n_atoms = len(ts)
    masses = ts.get_masses()

    if not np.all(masses):
        raise ValueError('Zero mass encountered in one or more of '
                         'the vibrated atoms. Use Atoms.set_masses()'
                         ' to set all masses to non-zero values.')
    mass_weights = np.repeat(masses**-0.5, 3)

    positions = ts.get_positions() - ts.get_center_of_mass()
    _, vectors_inertia = ts.get_moments_of_inertia(vectors=True)
    vectors_transrot = np.zeros((6, n_atoms, 3))
    vectors_transrot[0, :, 0] = 1
    vectors_transrot[1, :, 1] = 1
    vectors_transrot[2, :, 2] = 1
    vectors_transrot[3] = positions @ vectors_inertia[[1]].T @ vectors_inertia[[2]] - positions @ vectors_inertia[[2]].T @ vectors_inertia[[1]]
    vectors_transrot[4] = positions @ vectors_inertia[[2]].T @ vectors_inertia[[0]] - positions @ vectors_inertia[[0]].T @ vectors_inertia[[2]]
    vectors_transrot[5] = positions @ vectors_inertia[[0]].T @ vectors_inertia[[1]] - positions @ vectors_inertia[[1]].T @ vectors_inertia[[0]]
    vectors_transrot = vectors_transrot.reshape((6, n_atoms * 3))
    vectors_transrot = vectors_transrot / mass_weights
    vectors_transrot, _ = np.linalg.qr(vectors_transrot.T)
    vectors_transrot = vectors_transrot.T
    proj = np.eye(n_atoms * 3) - vectors_transrot.T @ vectors_transrot

    omega2, vectors = np.linalg.eigh(
        proj.T @ (
            mass_weights
            * get_hessian_2d(ts)
            * mass_weights[:, np.newaxis])
        @ proj)

    unit_conversion = units._hbar * units.m / sqrt(units._e * units._amu)
    energies = unit_conversion * omega2.astype(complex)**0.5

    modes = vectors.T.reshape(n_atoms * 3, n_atoms, 3)
    modes = modes * masses[np.newaxis, :, np.newaxis]**-0.5

    return (energies, modes)


def setup_analysis(
                   init_path,
                   ts,
                   highest_geodesic_ts_traj_file,
                   highest_geodesic_geom_file='highest_geodesic.xyz',
                  ):
    ts = read(xyz_ts)
    ts.calc = calc()
    sella_wrapper(ts, traj_file, sella_order=1)

    ts_e = ts.get_potential_energy()
    r_e = reactant.get_potential_energy()
    p_e = product.get_potential_energy()

    energy_list = [r_e]
    for image in images[1:-1]:
        image.calc = calc()
        energy_list.append(image.get_potential_energy())
    energy_list.append(p_e)
    write(highest_geodesic_geom_file, images[np.argmax(energy_list)])

    highest_geodesic = copy.deepcopy(images[np.argmax(energy_list)])
    highest_geodesic.calc = calc()
    sella_wrapper(highest_geodesic, highest_geodesic_ts_traj_file, sella_order=1)

    eigs, _ = _energies_and_modes(ts)
    print('\n', eigs)

    eigs, _ = _energies_and_modes(highest_geodesic)
    print('\n', eigs)

    print(f'Forward_barrier: {ts_e - r_e}')
    print(f'Reverse_barrier: {ts_e - p_e}')
