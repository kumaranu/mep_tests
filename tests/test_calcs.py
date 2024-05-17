import pytest
from calcs import load_config, calc_mace, calc
from mace.calculators import MACECalculator
from newtonnet.utils.ase_interface import MLAseCalculator


def test_load_config():
    config = load_config('tests/config.yml')

    assert isinstance(config, dict)

    assert 'mace' in config or 'newtonnet' in config, "Config should have either 'mace' or 'newtonnet' key"

    if 'mace' in config:
        assert 'model_path' in config['mace'], "'mace' should have a 'model_path' key"

    if 'newtonnet' in config:
        assert 'model_path' in config['newtonnet'], "'newtonnet' should have a 'model_path' key"
        assert 'config_path' in config['newtonnet'], "'newtonnet' should have a 'config_path' key"


def test_calc_mace():
    ml_calculator = calc_mace()
    assert isinstance(ml_calculator, MACECalculator)


def test_calc():
    ml_calculator = calc()
    assert isinstance(ml_calculator, MLAseCalculator)
