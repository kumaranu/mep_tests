import os
from calcs import calc
from ase.io import read, write
from typing import Any, Dict, Iterator, List, Sequence, Tuple, TypeVar, Union
from sella_wrapper import sella_wrapper


def setup_images(
                 logdir=None,
                 xyz_r_p=None,
                 xyz_ts=None,
                 N_intermediate=40,
                ):
    reactant = read(xyz_r_p, index='0')
    reactant.calc = calc()
    sella_wrapper(reactant, traj_file=logdir + '/reactant_opt.traj', sella_order=0)

    product = read(xyz_r_p, index='1')
    product.calc = calc()
    sella_wrapper(product, traj_file=logdir + '/product_opt.traj', sella_order=0)

    os.makedirs(logdir, exist_ok=True)
    os.chdir(logdir)
    write("r_p.xyz", [reactant.copy(), product.copy()])
    os.system('geodesic_interpolate r_p.xyz --output output.xyz --nimages ' + str(N_intermediate))
    images = read('output.xyz', index=':')

    for image in images:
        image.calc = calc()
        mlcalculator = calc()
        mlcalculator.calculate(image)
 
        energy = mlcalculator.results['energy']
        forces = mlcalculator.results['forces']
 
        image.info['energy'] = energy
        #image.arrays['forces'] = forces
 
    write('geodesic_path.xyz', images)

    return images
