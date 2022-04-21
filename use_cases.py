""" This file contains use cases. It can be executed as script and prints results to the console."""

from gpcrml_lib.data import ReceptorLibrary, Descriptors
from gpcrml_lib.machine_learning import DTCPredictions

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
print('\nPDBs with phi angles:\n\n', phi[:10])
print('\nPDBs with psi angles:\n\n', psi[:10])

# predict state of PDB files
print()
print('##########################')
print('# Calculate predictions. #')
print('##########################\n')
predictions = DTCPredictions(phi, psi)
predictions.predict()

# get phi and psi angles for an MD simulation
print()
print('################################################')
print('# Get phi and psi angles for an MD simulation. #')
print('################################################\n')
test_data_directory = 'gpcrml_lib/test_data'
rl = ReceptorLibrary()
rl.add_receptor(topology_path=test_data_directory + '/P29274.pdb', uniprotid='P29274',
                trajectory_path=test_data_directory + '/P29274.dcd')
d = Descriptors(rl)
phi = d.get_dihedrals('phi')
psi = d.get_dihedrals('psi')
print('\nMolecular dynamics simulation frames with phi angles:\n\n', phi[:10])
print('\nMolecular dynamics simulation frames with psi angles:\n\n', psi[:10])

# calculate predictions of MD simulation
print()
print('##########################')
print('# Calculate predictions. #')
print('##########################\n')
predictions = DTC(phi, psi)
predictions.predict()
print('Plotting dihedrals...')
predictions.plot_dihedrals()
print('Plotting prediction results...')
predictions.plot_predictions()
