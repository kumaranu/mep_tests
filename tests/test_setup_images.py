import os
import pytest
import shutil
from ase import Atoms
from ase.io import write
from setup_images import setup_images


@pytest.fixture
def setup_test_environment(tmp_path):
    # Create temporary directory
    logdir = tmp_path / "log"
    os.makedirs(logdir, exist_ok=True)

    # Create a mock XYZ file with reactant and product structures
    xyz_r_p = tmp_path / "r_p.xyz"

    reactant = Atoms(
        symbols='CCHHCHH',
        positions=[
            [1.4835950817281542, -1.0145410211301968, -0.13209027203235943],
            [0.8409564131524673, 0.018549610257914483, -0.07338809662321308],
            [-0.6399757891931867, 0.01763740851518944, 0.0581573443268891],
            [-1.0005576455546672, 1.0430257532387608, 0.22197240310602892],
            [1.402180736662139, 0.944112416574632, -0.12179540364365492],
            [-1.1216961389434357, -0.3883639833876232, -0.8769102842015071],
            [-0.9645026578514683, -0.6204201840686793, 0.9240543090678239]
        ]
    )

    product = Atoms(
        symbols='CCHHCHH',
        positions=[
            [1.348003553501624, 0.4819311116778978, 0.2752537177143993],
            [0.2386618286631742, -0.3433222966734429, 0.37705518940917926],
            [-0.9741307940518336, 0.07686022294949588, 0.08710778043683955],
            [-1.8314843503320921, -0.5547344604780035, 0.1639037492534953],
            [0.3801391040059668, -1.3793340533058087, 0.71035902765307],
            [1.9296265384257907, 0.622088341468767, 1.0901733942191298],
            [-1.090815880212625, 1.0965111343610956, -0.23791518420660265]
        ]
    )

    write(xyz_r_p, [reactant, product])

    return logdir, xyz_r_p


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
