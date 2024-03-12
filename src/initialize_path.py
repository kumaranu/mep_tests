import os
from ase.io import read
from calcs import calc
from ase.optimize import BFGS, ODE12r


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
