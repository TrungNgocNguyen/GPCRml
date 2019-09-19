from gpcrml_lib.data import ReceptorLibrary, Descriptors

library_directory = '/home/david/Documents/pycharm_projects/GPCRml/gpcrml_lib/test_data'

rl = ReceptorLibrary()
rl.add_receptors_from_pdb(pdb_codes=['3rze', '3pbl', '4n6h'], library_directory=library_directory)
rl.add_receptor(topology_path=library_directory + '/P29274.pdb', uniprotid='P29274', trajectory_path=library_directory +
                '/P29274.dcd')
d = Descriptors(rl)
phi = d.get_dihedrals('phi')
psi = d.get_dihedrals('psi')
