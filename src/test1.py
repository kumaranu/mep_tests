

from ase.optimize import BFGS
from ase.calculators.emt import EMT
from ase.mep.neb import NEBOptimizer
from neb_wrapper import run_neb_method

if __name__ == "__main__":
    # Define multiple sets of inputs
    input_sets = [
        #('aseneb', BFGS, None, None),
        #('improvedtangent', BFGS, None, None),
        #('spline', BFGS, None, None),
        #('string', BFGS, None, None),
        ('aseneb', NEBOptimizer, None, None),
        ('improvedtangent', NEBOptimizer, None, None),
        ('spline', NEBOptimizer, None, None),
        ('string', NEBOptimizer, None, None),
    ]

    # Set the test directory (you can modify this based on your needs)
    testdir = '/global/home/users/kumaranu/Documents/gpu_jobs'

    # Iterate over input sets and call run_neb_method
    for input_set in input_sets:
        run_neb_method(
            method=input_set[0],
            optimizer=input_set[1],
            precon=input_set[2],
            optmethod=input_set[3],
            testdir=testdir,
        )
