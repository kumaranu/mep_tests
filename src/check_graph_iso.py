import glob, os
from ase.io import read

ref_dir = '/global/cfs/cdirs/m2834/kumaranu/neb_nn_inputs'

for i in range(265):
    logdir = os.path.join(ref_dir, f'{i:03}')
    if os.path.exists(logdir):
        reactant_traj_file = os.path.join(logdir, 'reactant_opt.traj')
        if os.path.exists(reactant_traj_file):
            traj_atoms = read(reactant_traj_file, index=':', format="traj")
        else:
            print(f"Reactant's traj file could not be found for index {i:03}")
        #for atoms in traj_atoms:
        #    print(i, atoms)
    else:
        print(f'Path could not be found for index {i:03}.')

