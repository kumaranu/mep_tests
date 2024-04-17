import os
import torch
from mace.calculators import MACECalculator
from newtonnet.utils.ase_interface import MLAseCalculator
from ase import Atoms


def calc_mace():

    mlcalculator = MACECalculator(
        model_paths='/global/home/users/kumaranu/Documents/gpu_jobs/MACE_model.model',
        #model_paths='/global/home/users/kumaranu/Documents/mep_tests/src/MACE_model_cpu.model',
        device='cuda',
        #device='cpu',
        default_dtype="float64",
    )
    return mlcalculator


def calc():
    use_cuda = torch.cuda.is_available()
    device = 'cuda' if use_cuda else 'cpu'
    mlcalculator = MLAseCalculator(
        model_path=[
                    os.path.join(
                        os.path.expanduser("~"),
                        'Documents/NewtonNet/example/predict/training_52/models/best_model_state.tar'
                    ),
                    #'/global/home/users/kumaranu/Documents/NewtonNet/example/predict/training_53/models/best_model_state.tar',
                    #'/global/home/users/kumaranu/Documents/NewtonNet/example/predict/training_54/models/best_model_state.tar',
                    #'/global/home/users/kumaranu/Documents/NewtonNet/example/predict/training_55/models/best_model_state.tar',
                   ],    # path to model file, str or list of str
        settings_path=[
                       os.path.join(
                           os.path.expanduser("~"),
                           'Documents/NewtonNet/example/predict/training_52/run_scripts/config0.yml'
                       ),
                       #'/global/home/users/kumaranu/Documents/NewtonNet/example/predict/training_53/run_scripts/config2.yml',
                       #'/global/home/users/kumaranu/Documents/NewtonNet/example/predict/training_54/run_scripts/config1.yml',
                       #'/global/home/users/kumaranu/Documents/NewtonNet/example/predict/training_55/run_scripts/config3.yml',
                      ],    # path to configuration file, str or list of str
        hess_method=None,    # method to calculate hessians. 'autograd', 'fwd_diff', 'cnt_diff', or None (default: 'autograd')
        # hess_precision=1e-5,    # hessian gradient calculation precision for 'fwd_diff' and 'cnt_diff', ignored otherwise (default: None)
        disagreement='std',    # method to calculate disagreement among models. 'std', 'std_outlierremoval', 'range':, 'values', or None (default: 'std')
        device=device
    )
    return mlcalculator

'''
h2 = Atoms(numbers=[1, 1], positions=[[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
mlcalculator.calculate(h2)

print(mlcalculator.results['energy'])    # mean of calculated molecular energies, shape (1,)
print(mlcalculator.results['forces'])    # mean of calculated atomic forces, shape (n_atom, 3)
print(mlcalculator.results['hessian'])    # mean of calculated atomic Hessian, shape (n_atom, 3, n_atom, 3)
print(mlcalculator.results['energy_disagreement'])    # disagreement among calculated molecular 
'''

