import sys
import ase
import toml
from neb_wrapper import run_neb_method


def main(config_file: str) -> None:
    """
    Run NEB method using the provided configuration file.

    Args:
        config_file (str): Path to the configuration file.
    """
    with open(config_file, 'r') as f:
        inputs = toml.load(f)
    
    input_sets = inputs['input_sets']

    for input_set in input_sets:
        optimizer_name = input_set['optimizer']
        # ase_optimizer = getattr(ase.mep.neb, optimizer_name)
        ase_optimizer = getattr(ase.neb, optimizer_name)

        run_neb_method(
            method=input_set['method'],
            optimizer=ase_optimizer,
            opt_method=None,
            precon=None,
            logdir=inputs['input_paths']['logdir'],
            xyz_r_p=inputs['input_paths']['xyz_r_p'],
            xyz_ts=inputs['input_paths']['xyz_ts'],
            n_intermediate=inputs['geodesic_inputs']['N_intermediate'],
        )


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py config_file.toml")
        sys.exit(1)

    config_file = sys.argv[1]
    main(config_file)
