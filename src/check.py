import glob, os, sys

logdir = '/global/scratch/users/kumaranu/saved_files/neb_nn_inputs'

os.chdir(logdir)
for i in range(265):
    if os.path.exists(f'{i:03}'):
        os.chdir(f'{i:03}')
        if os.path.exists('neb_band_aseneb_NEBOptimizer_None.txt'):
            INF = open('neb_band_aseneb_NEBOptimizer_None.txt', 'r').readlines()
            #if int(INF[-1].split()[1]) > 200:
            print(f'{i:03}: {INF[-1].split()[1]}')
        else:
            print(f'{i:03}: No NEB file.')
        os.chdir('../')
    else:
        print(f'{i:03} does not exist.')

