""" This file contains use cases. It can be executed as script and prints results to the console."""

from gpcrml_lib.data import ReceptorLibrary, Descriptors

# get phi and psi angles for a set of PDB codes
print()
print('##################################################')
print('# Get phi and psi angles for a set of PDB codes. #')
print('##################################################\n')
test_data_directory = 'gpcrml_lib/test_data'
rl = ReceptorLibrary()
rl.add_receptors_from_pdb(pdb_codes=['3rze', '3pbl', '4n6h'], library_directory=test_data_directory)
d = Descriptors(rl)
phi = d.get_dihedrals('phi')
psi = d.get_dihedrals('psi')
print('\nFirst 10 residues with phi angles:\n\n', phi[:10])
print('\nFirst 10 residues with psi angles:\n\n', psi[:10])

# get phi and psi angles for an MD simulation
print()
print('################################################\n')
print('# Get phi and psi angles for an MD simulation. #\n')
print('################################################\n')
test_data_directory = 'gpcrml_lib/test_data'
rl = ReceptorLibrary()
rl.add_receptor(topology_path=test_data_directory + '/P29274.pdb', uniprotid='P29274',
                trajectory_path=test_data_directory + '/P29274.dcd')
d = Descriptors(rl)
phi = d.get_dihedrals('phi')
psi = d.get_dihedrals('psi')
print('\nFirst 10 residues with phi angles:\n\n', phi[:10])
print('\nFirst 10 residues with psi angles:\n\n', psi[:10])
