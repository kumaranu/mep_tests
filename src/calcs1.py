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
    mlcalculator = MLAseCalculator(
        model_path=[
                    '/global/home/users/kumaranu/Documents/NewtonNet/example/predict/training_52/models/best_model_state.tar',
                   ],    # path to model file, str or list of str
        settings_path=[
                       '/global/home/users/kumaranu/Documents/NewtonNet/example/predict/training_52/run_scripts/config0.yml',
                      ],    # path to configuration file, str or list of str
        hess_method=None,    # method to calculate hessians. 'autograd', 'fwd_diff', 'cnt_diff', or None (default: 'autograd')
        # hess_precision=1e-5,    # hessian gradient calculation precision for 'fwd_diff' and 'cnt_diff', ignored otherwise (default: None)
        disagreement='std',    # method to calculate disagreement among models. 'std', 'std_outlierremoval', 'range':, 'values', or None (default: 'std')
        #device='cuda'   # 'cpu' or list of cuda
        device='cpu'   # 'cpu' or list of cuda
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

