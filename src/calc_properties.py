from ase.io import read, write
from calcs import calc

# Read molecular configurations from XYZ file
atoms_list = read('/global/scratch/users/kumaranu/saved_files/neb_nn_inputs/263/output.xyz', index=':')

# Create a calculator instance (e.g., EMT)
mlcalculator = calc()

# Iterate over each molecular configuration
for atoms in atoms_list:
    mlcalculator.calculate(atoms)
    
    # Calculate energies and forces
    energy = mlcalculator.results['energy']
    forces = mlcalculator.results['forces']
    
    # Store energies and forces in the atoms object
    atoms.info['energy'] = energy
    atoms.arrays['forces'] = forces

# Write the updated atoms object with energies and forces to a new XYZ file
write('/global/scratch/users/kumaranu/saved_files/neb_nn_inputs/263/output1.xyz', atoms_list)


'''
h2 = Atoms(numbers=[1, 1], positions=[[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
mlcalculator.calculate(h2)

print(mlcalculator.results['energy'])    # mean of calculated molecular energies, shape (1,)
print(mlcalculator.results['forces'])    # mean of calculated atomic forces, shape (n_atom, 3)
print(mlcalculator.results['hessian'])    # mean of calculated atomic Hessian, shape (n_atom, 3, n_atom, 3)
print(mlcalculator.results['energy_disagreement'])    # disagreement among calculated molecular 
'''



