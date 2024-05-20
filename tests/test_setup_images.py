import os
import pytest
from ase import Atoms
from ase.io import write
from setup_images import setup_images


def test_setup_images(setup_test_environment):
    logdir, xyz_r_p = setup_test_environment

    # Call the setup_images function
    images = setup_images(logdir=str(logdir), xyz_r_p=str(xyz_r_p), n_intermediate=2)

    # Check that images were returned
    assert len(images) > 0, "No images were generated"

    # Verify output files were created
    assert os.path.isfile(logdir / 'reactant_opt.traj'), "Reactant optimization file not found"
    assert os.path.isfile(logdir / 'product_opt.traj'), "Product optimization file not found"
    assert os.path.isfile(logdir / 'r_p.xyz'), "Reactant-Product file not found"
    assert os.path.isfile(logdir / 'output.xyz'), "Intermediate images file not found"
    assert os.path.isfile(logdir / 'geodesic_path.xyz'), "Geodesic path file not found"

    # Check energies and forces
    for image in images:
        assert 'energy' in image.info, "Energy not found in image info"
        assert 'forces' in image.arrays, "Forces not found in image arrays"


if __name__ == "__main__":
    pytest.main()
