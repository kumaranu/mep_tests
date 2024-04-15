
import sys
import ase
import toml
from neb_wrapper import run_neb_method


if __name__ == "__main__":
    with open(sys.argv[1], 'r') as f:
        inputs = toml.load(f)
    
    input_sets = inputs['input_sets']
    logdir = inputs['input_paths']['logdir']

    for input_set in input_sets:
        optimizer_name = input_set['optimizer']
        # ase_optimizer = getattr(ase.mep.neb, optimizer_name)
        ase_optimizer = getattr(ase.neb, optimizer_name)

        run_neb_method(
            method=input_set['method'],
            optimizer=ase_optimizer,
            precon=None,
            optmethod=None,
            logdir=inputs['input_paths']['logdir'],
            xyz_r_p=inputs['input_paths']['xyz_r_p'],
            xyz_ts=inputs['input_paths']['xyz_ts'],
            N_intermediate=inputs['geodesic_inputs']['N_intermediate'],
        )

