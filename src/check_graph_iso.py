import glob
import os
from ase.io import read
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN
from networkx.algorithms.graph_hashing import weisfeiler_lehman_graph_hash


def add_specie_suffix(graph):
    """
    Add a suffix to each node's 'specie' attribute in the graph.

    Args:
        graph (networkx.Graph): The graph.
    """
    for idx in graph.nodes():
        graph.nodes()[idx]["specie"] = graph.nodes()[idx]["specie"] + str(idx)


def get_graph_hash(graph):
    """
    Get the hash of the graph using the Weisfeiler-Lehman algorithm.

    Args:
        graph (networkx.Graph): The graph.

    Returns:
        str: The graph hash.
    """
    return weisfeiler_lehman_graph_hash(graph, node_attr='specie')


def create_molecule_graph(molecule):
    """
    Create a molecule graph using the OpenBabelNN strategy.

    Args:
        molecule (pymatgen.core.structure.Molecule): The molecule.

    Returns:
        pymatgen.analysis.graphs.MoleculeGraph: The molecule graph.
    """
    return MoleculeGraph.with_local_env_strategy(molecule, OpenBabelNN())


def compare_mols(molecule_1, molecule_2) -> bool:
    """
    Compare two molecules based on their graph structure.

    Args:
        molecule_1 (pymatgen.core.structure.Molecule): First molecule.
        molecule_2 (pymatgen.core.structure.Molecule): Second molecule.

    Returns:
        bool: True if the molecules have the same graph structure, False otherwise.
    """
    molgraph_1 = create_molecule_graph(molecule_1)
    molgraph_2 = create_molecule_graph(molecule_2)

    graph_1 = molgraph_1.graph.to_undirected()
    graph_2 = molgraph_2.graph.to_undirected()

    add_specie_suffix(graph_1)
    add_specie_suffix(graph_2)

    graph_1_hash = get_graph_hash(graph_1)
    graph_2_hash = get_graph_hash(graph_2)

    return graph_1_hash == graph_2_hash


def check_graph_iso(ref_dir):
    num_directories = len(glob.glob(ref_dir + '/???/'))
    for i in range(num_directories):
        logdir = os.path.join(ref_dir, f'{i:03}')
        if os.path.exists(logdir):
            reactant_traj_file = os.path.join(logdir, 'reactant_opt.traj')
            if os.path.exists(reactant_traj_file):
                reactant_traj_atoms = read(reactant_traj_file, index=':', format="traj")
            else:
                print(f"Reactant's traj file could not be found for index {i:03}")
                continue
            product_traj_file = os.path.join(logdir, 'product_opt.traj')
            if os.path.exists(product_traj_file):
                pdt_traj_atoms = read(product_traj_file, index=':', format="traj")
            else:
                print(f"Product's traj file could not be found for index {i:03}")
                continue
            traj_pymat_mol1 = AseAtomsAdaptor.get_molecule(reactant_traj_atoms[0])
            traj_pymat_mol2 = AseAtomsAdaptor.get_molecule(reactant_traj_atoms[-1])
            bool_reactant = compare_mols(traj_pymat_mol1, traj_pymat_mol2)
            traj_pymat_mol1 = AseAtomsAdaptor.get_molecule(pdt_traj_atoms[0])
            traj_pymat_mol2 = AseAtomsAdaptor.get_molecule(pdt_traj_atoms[-1])
            bool_product = compare_mols(traj_pymat_mol1, traj_pymat_mol2)
            if bool_reactant and bool_product:
                continue
            else:
                print(f'bool_reactant: {bool_reactant}, bool_product: {bool_product} for index {i:03}')
        else:
            print(f'Path could not be found for index {i:03}.')


if __name__ == '__main__':
    # ref_dir = '/global/cfs/cdirs/m2834/kumaranu/neb_nn_inputs'
    ref_dir = '/global/cfs/cdirs/m2834/kumaranu/neb_nn_inputs'
    # ref_dir = '/home/kumaranu/Downloads/neb_nn_inputs'
    check_graph_iso(ref_dir)
