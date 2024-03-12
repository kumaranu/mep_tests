
import os
import matplotlib.pyplot as plt
from initialize_path import setup_images
from ase.neb import NEB
from ase.neb import NEBTools


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
