import os
import pytest
from ase import Atoms
from ase.io import read, write
from setup_images import setup_images
from setup_images import geodesic_interpolate_wrapper


def test_geodesic_interpolate_wrapper(setup_test_environment):
    logdir, xyz_r_p = setup_test_environment
    atoms_object = read(xyz_r_p, index=':')
    symbols, smoother_path = geodesic_interpolate_wrapper(
        atoms_object
    )
    # assert output == 1
    # assert symbols == 1
    assert smoother_path[1][0][0] == pytest.approx(
        1.36055556030,
        abs=1e-3,
    )


def test_setup_images(setup_test_environment):
    logdir, xyz_r_p = setup_test_environment

    # Call the setup_images function
    images = setup_images(
        logdir=str(logdir),
        xyz_r_p=str(xyz_r_p),
        n_intermediate=2
    )

    # Check that images were returned
    assert len(images) > 0, "No images were generated"

    # Verify output files were created
    assert os.path.isfile(logdir / 'reactant_opt.traj'), "Reactant optimization file not found"
    assert os.path.isfile(logdir / 'product_opt.traj'), "Product optimization file not found"
    assert os.path.isfile(logdir / 'r_p.xyz'), "Reactant-Product file not found"

    assert images[1].get_positions()[0][0] == pytest.approx(
        1.439202,
        abs=1e-2,
    )

    # Check energies and forces
    for image in images:
        assert 'energy' in image.info, "Energy not found in image info"
        assert 'forces' in image.arrays, "Forces not found in image arrays"

    assert images[1].get_potential_energy() == pytest.approx(
        -29.9266674,
        abs=1,
    )
    "Error in first intermediate image's energy for the geodesic path"

    assert images[0].get_potential_energy() == pytest.approx(
        -32.98172029,
        abs=1,
    )
    "Error in reactant energy prediction for geodesic path."

    assert images[-1].get_potential_energy() == pytest.approx(
        -25.1172894,
        abs=1,
    )
    "Error in product energy prediction for geodesic path."


if __name__ == "__main__":
    pytest.main()
