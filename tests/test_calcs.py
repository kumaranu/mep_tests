import pytest
from ase.io import read
from calcs import load_config, calc_mace, calc
from mace.calculators import MACECalculator
from newtonnet.utils.ase_interface import MLAseCalculator
#from conftest import setup_test_environment


def test_load_config():
    config = load_config('tests/config.yml')

    assert isinstance(config, dict)

    assert 'mace' in config or 'newtonnet' in config, "Config should have either 'mace' or 'newtonnet' key"

    if 'mace' in config:
        assert 'model_path' in config['mace'], "'mace' should have a 'model_path' key"

    if 'newtonnet' in config:
        assert 'model_path' in config['newtonnet'], "'newtonnet' should have a 'model_path' key"
        assert 'config_path' in config['newtonnet'], "'newtonnet' should have a 'config_path' key"


def test_calc_mace(setup_test_environment):
    ml_calculator = calc_mace()
    assert isinstance(ml_calculator, MACECalculator)


def test_calc(setup_test_environment):
    ml_calculator = calc()
    assert isinstance(ml_calculator, MLAseCalculator)

    logdir, xyz_r_p = setup_test_environment
    print('xyz_r_p:', xyz_r_p)
    atoms_object = read(xyz_r_p)
    ml_calculator.calculate(atoms_object)

    assert ml_calculator.results['energy'] == pytest.approx(
        -20.7079663,
        abs=1e-5,
    )
    "error in energy"

    assert ml_calculator.results['forces'][0, 0] == pytest.approx(
        5.6846989,
        abs=1e-5,
    )
    "error in forces"

    '''
    !!!!! HESSIAN TESTS FAIL RIGHT NOW !!!!!
    !!!!! NEED TO TALK TO Eric !!!!!
    assert sum(sum(sum(ml_calculator.results['hessian']))) == pytest.approx(
        5.6846989,
        abs=1e-5,
    )
    "error in Hessian"
    '''
