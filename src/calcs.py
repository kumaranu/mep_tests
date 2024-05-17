import os
import yaml
import torch
from mace.calculators import MACECalculator
from newtonnet.utils.ase_interface import MLAseCalculator
from typing import Any, Dict


def load_config(config_path: str = 'tests/config.yml') -> Dict[str, Any]:
    """
    Load the configuration file.

    Parameters
    ----------
    config_path : str, optional
        Path to the configuration file, by default 'config.yml'

    Returns
    -------
    Dict[str, Any]
        Configuration dictionary.
    """
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)
    return config


def calc_mace() -> MACECalculator:
    """
    Initialize and return a MACECalculator instance.

    Returns
    -------
    MACECalculator
        An instance of the MACECalculator initialized with the specified model path and device.
    """
    config = load_config()
    model_path = config['mace']['model_path']
    device = 'cuda' if torch.cuda.is_available() else 'cpu'

    ml_calculator = MACECalculator(
        model_paths=model_path,
        device=device,
        default_dtype="float64",
    )
    return ml_calculator


def calc() -> MLAseCalculator:
    """
    Initialize and return an MLAseCalculator instance.

    Returns
    -------
    MLAseCalculator
        An instance of the MLAseCalculator initialized with the specified model path, config path, and device.
    """
    config = load_config()
    model_path = config['newtonnet']['model_path']
    config_path = config['newtonnet']['config_path']
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    ml_calculator = MLAseCalculator(
        model_path=[model_path],
        settings_path=[config_path],
        disagreement='std',
        device=device
    )
    return ml_calculator
