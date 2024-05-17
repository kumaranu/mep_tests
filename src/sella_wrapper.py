from ase.io import write
from ase.io import Trajectory
from sella import Sella
from typing import Any, Dict, Iterator, List, Sequence, Tuple, TypeVar, Union


def sella_wrapper(
                  atoms_object,
                  traj_file=None,
                  sella_order=0,
                  use_internal=True,
                  traj_log_interval=2,
                  fmax_cutoff=1e-3,
                  max_steps=1000
                 ):
    if traj_file:
        traj = Trajectory(
                          traj_file,
                          'w',
                          atoms_object,
                         )
    qn = Sella(
               atoms_object,
               order=sella_order,
               internal=use_internal,
              )
    if traj_file:
        qn.attach(
                  traj.write,
                  interval=traj_log_interval,
                 )
    qn.run(
           fmax=fmax_cutoff,
           steps=max_steps,
          )
    if traj_file:
        traj.close()
