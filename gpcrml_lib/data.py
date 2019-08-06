from urllib.request import urlopen
import json
from collections import OrderedDict
from Bio import pairwise2
from checker import special_character, len_is_four
from timeit import default_timer as timer


class GetReceptorMapping:
    """A descriptor extracted from a PDB file or MD trajectory, mapped to BWN numbering scheme"""

    def __init__(self, topology_dir, trajectory_dir, identifier):
        self.topology_dir = topology_dir
        self.trajectory_dir = trajectory_dir
        self.identifier = identifier

        self.protein_name = None
        self.pdb_code = None
        self.pdb_marker = None
        self.preferred_chain = None

        self.sequence_dict = None
        self.pdb_sequence_dict = None
        self.pdb_generic_numbers_dict = None
        self.pdb_generic_numbers_dict_marker = None

    def align_with_generic_numbers(self):
        # Helper functions
        # Get residue and bwn list based on protein_name (uniprotid name)
        def get_protein_name(identifier):
            """
            Checker if identifier is 4 characters long, if yes get protein_name, pdb code, pref_chain with help of
            gpcrdb, else it should be the protein_name
            :param identifier: self.identifier
            :return: write protein_name in class attribute
            """
            if len_is_four(identifier) and not special_character(identifier):
                print('Identifier is PDB code "{}"...'.format(identifier))
                response = urlopen('https://gpcrdb.org/services/structure/{}'.format(identifier))
                protein = json.loads(response.read().decode('utf-8'))
                self.preferred_chain = protein.get('preferred_chain', '')
                self.pdb_code = identifier
                self.protein_name = protein.get('protein', '')
                self.pdb = True
            else:
                print('Identifier is protein name "{}"...'.format(identifier))
                self.pdb = False
                if special_character(identifier, '_'):
                    self.protein_name = identifier
                else:
                    print("The specified identifier '{}' needs to be a PDB code (e.g. '3V2Y') ".format(identifier),
                          "or a protein name (e.g. 's1pr1_human')")

        def get_sequence_dict(protein_name):
            """"""
            print('Retrieving sequence for "{}"...'.format(protein_name))
            response = urlopen('https://gpcrdb.org/services/residues/{}'.format(protein_name))
            protein = json.loads(response.read().decode('utf-8'))
            dict = OrderedDict()
            for residue in protein:
                dict[residue['sequence_number']] = [residue['amino_acid'], residue['display_generic_number']]
            self.sequence_dict = dict

        def get_pdb_sequence_dict(preferred_chain, topology_dir):
            amino_acid_dict = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E',
                               'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F',
                               'PRO': 'P', 'SEC': 'U', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}
            pdb_sequence_dict = OrderedDict()
            print('Creating sequence dict of "{}" ...'.format(topology_dir))
            with open(topology_dir, 'r') as file:
                for line in file.readlines():
                    if 'ATOM' in line[:6]:
                        if self.pdb_marker:
                            if line[21:22] == preferred_chain:
                                if line[17:20] in amino_acid_dict.keys():
                                    if line[13:16].strip() == 'N':
                                        pdb_sequence_dict[int(line[22:26].strip())] = amino_acid_dict[line[17:20]]
                        else:
                            if line[17:20] in amino_acid_dict.keys():
                                if line[13:16].strip() == 'N':
                                    pdb_sequence_dict[int(line[22:26].strip())] = amino_acid_dict[line[17:20]]
            self.pdb_sequence_dict = pdb_sequence_dict

        def assign_generic_numbers_to_pdb(sequence_dict, pdb_sequence_dict):
            print('Assigning generic numbers.')
            sequence = ''.join([sequence_dict[residue][0] for residue in sequence_dict.keys()])
            pdb_sequence = ''.join([pdb_sequence_dict[residue] for residue in pdb_sequence_dict.keys()])
            alignment = pairwise2.align.globalxs(sequence, pdb_sequence, -10, 0)[0]
            pdb_generic_numbers_dict = {}
            pdb_sequence_ids = list(pdb_sequence_dict.keys())
            counter = 1
            pdb_counter = 0
            for residue, pdb_residue in zip(alignment[0], alignment[1]):
                if residue != '-' and pdb_residue != '-':
                    pdb_generic_numbers_dict[pdb_sequence_ids[pdb_counter]] = sequence_dict[counter][1]
                    counter += 1
                    pdb_counter += 1
                else:
                    if residue != '-':
                        counter += 1
                    if pdb_residue != '-':
                        pdb_counter += 1
            self.pdb_generic_numbers_dict = pdb_generic_numbers_dict
            print('Process finished. Dictionary accessible under "your_instance.pdb_generic_numbers_dict" \n')


        get_protein_name(self.identifier)
        get_sequence_dict(self.protein_name)
        get_pdb_sequence_dict(self.preferred_chain, self.topology_dir)
        assign_generic_numbers_to_pdb(self.sequence_dict, self.pdb_sequence_dict)

    def get_dihedral(self):
        pass

    def get_descriptor_x(self):
        pass



# START MY TIMER
start = timer()

s1pr1 = GetReceptorMapping("/home/gocky/python/GPCR-ML_data/pdbs/3V2Y.pdb", "", "3V2Y")
s1pr1.align_with_generic_numbers()

s1pr2 = GetReceptorMapping("/mdspace/gocky/S1PR2_project/homology_models_e28_docked/s1pr2.pdb", "", "s1pr2_human")
s1pr2.align_with_generic_numbers()

# STOP MY TIMER
elapsed_time = timer() - start # in seconds
print(elapsed_time)


