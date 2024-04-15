import os
from calcs import calc
from ase.io import write
import matplotlib.pyplot as plt
from setup_images import setup_images
from ase.neb import NEB
from ase.neb import NEBTools


def run_neb_method(
                   method,
                   optimizer,
                   precon,
                   optmethod,
                   logdir=None,
                   xyz_r_p=None,
                   xyz_ts=None,
                   N_intermediate=20,
                  ):
    images = setup_images(logdir, xyz_r_p, xyz_ts, N_intermediate)
    
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
              climb=True,
              precon=precon,
              remove_rotation_and_translation=True,
              parallel=True,
             )

    os.makedirs(logdir, exist_ok=True)
    log_filename = f'neb_band_{method}_{optimizer.__name__}_{precon}.txt'

    logfile_path = os.path.join(logdir, log_filename)

    if optmethod is not None:
        opt = optimizer(mep, method=optmethod, logfile=logfile_path, verbose=2)
    else:
        opt = optimizer(mep, logfile=logfile_path, verbose=2)

    opt.attach(save_fmax_history, 1, mep)
    opt.run(fmax=1e-2, steps=400)
    # opt.run(fmax=1e-2, steps=1000)

    write(f'{logdir}/optimized_path_{method}_{optimizer.__name__}_{precon}.xyz', images)

    # nebtools = NEBTools(images)

    # Ef, dE = nebtools.get_barrier(fit=False)
    # print(f'{method},{optimizer.__name__},{precon} '
    #       f'=> Ef = {Ef:.3f}, dE = {dE:.3f}')

    # # Save the NEB band plot
    # output_filename = f'neb_band_{method}_{optimizer}_{precon}.png'
    # output_path = os.path.join(logdir, output_filename)

    # # Plot and save the NEB band
    # fig, ax = plt.subplots()
    # # nebtools.plot_band(ax=ax)

    # for image in initial_images_copy:
    #     image.calc = calc()
    #     ax.plot(image.get_potential_energy(), label='Initial Image')
    # for image in images:
    #     ax.plot(image.get_potential_energy(), label='Optimized Image')
    # ax.set_xlabel('Image Index')
    # ax.set_ylabel('Potential Energy (eV)')
    # ax.legend()

    # plt.savefig(output_path)
    # plt.close(fig)

