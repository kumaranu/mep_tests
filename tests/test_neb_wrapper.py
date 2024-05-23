import os

import numpy as np
import pytest
from ase.io import read
from neb_wrapper import run_neb_method
from ase.mep.neb import NEBOptimizer

def test_run_neb_method(tmp_path, setup_test_environment):
    logdir, xyz_r_p = setup_test_environment

    logdir = tmp_path / "logs"
    method = "aseneb"
    optimizer = NEBOptimizer
    opt_method = None
    precon = None
    n_intermediate = 10
    k = 0.1
    max_steps = 3
    fmax_cutoff = 1e-3

    images = run_neb_method(
        method=method,
        optimizer=optimizer,
        opt_method=opt_method,
        precon=precon,
        logdir=str(logdir),
        xyz_r_p=xyz_r_p,
        n_intermediate=n_intermediate,
        k=k,
        max_steps=max_steps,
        fmax_cutoff=fmax_cutoff,
    )
    print(f"neb_band_{method}_{optimizer.__name__}_{precon}.txt")

    assert images[0].positions[0][1] == pytest.approx(
        -0.85232544,
        abs=1e-2,
    )

    assert images[0].get_potential_energy() == pytest.approx(
        -32.9817202,
        abs=1e-2,
    ), "reactant potential energy"

    assert images[0].get_forces()[0, 1] == pytest.approx(
        -5.5567398e-05,
        abs=1e-5,
    ), "reactant potential forces"

    assert images[1].positions[0][1] == pytest.approx(
        -8.65177466e-01,
        abs=1e-2,
    )

    assert np.argmax(
        [
            image.get_potential_energy() for image in images
        ]
    ) == pytest.approx(
        4,
    ), "Index of the transition state"

    assert np.max(
        [
            image.get_potential_energy() for image in images
        ]
    ) == pytest.approx(
        -19.946616164,
        abs=1,
    ), "Potential energy of the transition state"

    assert images[
               np.argmax(
                   [
                       image.get_potential_energy() for image in images
                   ]
               )
           ].get_forces()[0, 1] == pytest.approx(
        -4.60381977,
        abs=1,
    ), "Force component in the transition state"

    assert images[-1].positions[0][1] == pytest.approx(
        -1.31163882,
        abs=1e-2,
    )

    assert os.path.exists(
        logdir / f"neb_band_{method}_{optimizer.__name__}_{precon}.txt"
    ), 'Could not find the optimization output file for NEB'

    assert os.path.exists(
        f'{logdir}/optimized_path_{method}_{optimizer.__name__}_{precon}.xyz',
    ), "Could not find the xyz file for converged NEB calculation."
