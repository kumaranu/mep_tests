import os
from ase.io import read
from ase.neb import NEB
from ase.neb import NEBTools
from ase.optimize import BFGS
import matplotlib.pyplot as plt
from ase.calculators.emt import EMT
from ase.mep.neb import NEBOptimizer
from ase.optimize import BFGS, ODE12r


def calc():
    from mace.calculators import MACECalculator
    mlcalculator = MACECalculator(
        model_paths='/global/home/users/kumaranu/Documents/gpu_jobs/MACE_model.model',
        device='cuda',
        default_dtype="float64",
    )
    return mlcalculator


def setup_images():
    N_intermediate = 40
    # xyz_r_p = '/global/home/users/kumaranu/trash/inputs/abcd/A_B.xyz'
    # xyz_r_p = '/global/home/users/kumaranu/trash/inputs/abcd/A_C+D.xyz'
    xyz_r_p = '/global/home/users/kumaranu/trash/inputs/abcd/B_C+D.xyz'
    reactant = read(xyz_r_p, index='0')
    reactant.calc = calc()
    product = read(xyz_r_p, index='1')
    product.calc = calc()
    os.makedirs('/tmp/geodesic_calc_dir', exist_ok=True)
    os.chdir('/tmp/geodesic_calc_dir')
    os.system('geodesic_interpolate ' + str(xyz_r_p) + ' --output output.xyz --nimages ' + str(N_intermediate))
    intermediate = read('output.xyz', index=':')

    qn = ODE12r(reactant)
    qn.run(fmax=1e-3, steps=1000)

    qn = ODE12r(product)
    qn.run(fmax=1e-3, steps=1000)

    images = [reactant]
    for image in intermediate:
        image.calc = calc()
        images.append(image)
    images.append(product)
    return images, None, None


def run_neb_method(method, optimizer, precon, optmethod, testdir='trash/'):
    images, _, _ = setup_images()

    fmax_history = []

    def save_fmax_history(mep):
        fmax_history.append(mep.get_residual())

    k = 0.1
    if precon == 'Exp':
        k = 0.01
    mep = NEB(
              images,
              k=k,
              method=method,
              climb=False,
              precon=precon,
              remove_rotation_and_translation=False,
              parallel=True,
             )

    os.makedirs(testdir, exist_ok=True)
    log_filename = f'neb_band_{method}_{optimizer.__name__}_{precon}.txt'
    logfile_path = os.path.join(testdir, log_filename)

    if optmethod is not None:
        opt = optimizer(mep, method=optmethod, logfile=logfile_path)
    else:
        opt = optimizer(mep, logfile=logfile_path)

    opt.attach(save_fmax_history, 1, mep)
    opt.run(fmax=1e-2, steps=1000)

    nebtools = NEBTools(images)

    Ef, dE = nebtools.get_barrier(fit=False)
    print(f'{method},{optimizer.__name__},{precon} '
          f'=> Ef = {Ef:.3f}, dE = {dE:.3f}')

    # Save the NEB band plot
    output_filename = f'neb_band_{method}_{optimizer.__name__}_{precon}.png'
    output_path = os.path.join(testdir, output_filename)

    # Plot and save the NEB band
    fig, ax = plt.subplots()
    nebtools.plot_band(ax=ax)
    plt.savefig(output_path)
    plt.close(fig)


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

