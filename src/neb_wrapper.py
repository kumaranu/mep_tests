import os
from typing import Optional, List, Union
from ase.io import write
from setup_images import setup_images
from ase.neb import NEB
from ase.optimize import Optimizer


def run_neb_method(
        method: str,
        optimizer: Optimizer,
        opt_method: Optional[str] = None,
        precon: Optional[str] = None,
        logdir: Optional[str] = None,
        xyz_r_p: Optional[str] = None,
        xyz_ts: Optional[str] = None,
        n_intermediate: Optional[int] = 20,
        k: Optional[float] = 0.1,
        max_steps: Optional[int] = 1000,
        fmax_cutoff: Optional[float] = 1e-2,
) -> None:
    """
    Run NEB method.

    Args:
        method (str): NEB method.
        optimizer: Optimizer function.
        precon (str, optional): Preconditioner method. Defaults to None.
        opt_method (str, Optimizer): Optimization method. Defaults to None.
        logdir (str, optional): Directory to save logs. Defaults to None.
        xyz_r_p (str, optional): Path to reactant and product XYZ files. Defaults to None.
        xyz_ts (str, optional): Path to transition state XYZ file. Defaults to None.
        n_intermediate (int, optional): Number of intermediate images. Defaults to 20.
        k (float, optional): force constant for the springs in NEB. Defaults to 0.1.
        max_steps (int, optional): maximum number of optimization steps allowed. Defaults to 1000.
        fmax_cutoff (float: optional): convergence cut-off criteria for the NEB optimization. Defaults to 1e-2.
    """
    images = setup_images(logdir, xyz_r_p, xyz_ts, n_intermediate)
    
    mep = NEB(
        images,
        k=k,
        method=method,
        climb=True,
        precon=precon,
        remove_rotation_and_translation=True,
        parallel=True,
    )

    os.makedirs(logdir, exist_ok=True)
    log_filename = f'neb_band_{method}_{optimizer.__name__}_{precon}.txt'

    logfile_path = os.path.join(logdir, log_filename)

    if opt_method is not None:
        opt = optimizer(mep, method=opt_method, logfile=logfile_path, verbose=2)
    else:
        opt = optimizer(mep, logfile=logfile_path, verbose=2)

    opt.run(fmax=fmax_cutoff, steps=max_steps)

    write(f'{logdir}/optimized_path_{method}_{optimizer.__name__}_{precon}.xyz', images)
