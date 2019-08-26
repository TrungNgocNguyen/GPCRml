import json
from urllib.request import urlopen
from collections import OrderedDict
from Bio import pairwise2
from checker import special_character, len_is_x, capital_letter_or_digit
from timeit import default_timer as timer

import MDAnalysis as mda
import pandas as pd

class GetReceptorMapping:
    """Map PDB sequence to BWN numbering scheme from GPCRdb"""

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

    def get_protein_name(self, identifier):
        """
        Checker if identifier is 4 characters long, if yes get protein_name, pdb code, pref_chain with help of
        gpcrdb, else it should be the protein_name
        :param identifier: self.identifier
        :return: write protein_name and other variables in instance attributes
        """
        if len_is_x(identifier, 4) and not special_character(identifier) and capital_letter_or_digit(identifier):
            print('Identifier is PDB code "{}"...'.format(identifier))
            response = urlopen('https://gpcrdb.org/services/structure/{}'.format(identifier))
            protein = json.loads(response.read().decode('utf-8'))
            self.preferred_chain = protein.get('preferred_chain', '')
            self.pdb_code = identifier
            self.protein_name = protein.get('protein', '')
            self.pdb_marker = True
        else:
            self.pdb_marker = False
            if special_character(identifier, '_'):
                self.protein_name = identifier
                print('Identifier is protein name "{}"...'.format(identifier))

            else:
                raise ValueError("The specified identifier '{}' needs to be a PDB code in capital letters (e.g. '3V2Y')"
                                 .format(identifier), " or a protein name in lowercase letters (e.g. 's1pr1_human')")

    def get_sequence_dict(self, protein_name):
        """
        Get sequence based on protein name from GPCRdb
        :param protein_name: retrieved based on PDB or through specified instance attribute (identifier)
        :return: Save dict in instance attribute "sequence_dict"
        """
        print('Retrieving sequence for "{}"...'.format(protein_name))
        response = urlopen('https://gpcrdb.org/services/residues/{}'.format(protein_name))
        protein = json.loads(response.read().decode('utf-8'))
        dictionary = OrderedDict()
        for residue in protein:
            dictionary[residue['sequence_number']] = [residue['amino_acid'], residue['display_generic_number']]
        self.sequence_dict = dictionary

    def get_pdb_sequence_dict(self, preferred_chain, topology_dir):
        replace = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E',
                   'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F',
                   'PRO': 'P', 'SEC': 'U', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
                   '<Residue': '', ' ': '', '>': ''}

        u = mda.Universe(topology_dir)

        if self.pdb_marker:
            u.select_atoms('segid '+preferred_chain)
        else:
            pass

        universe_list = list(u.atoms[:].residues)

        df = pd.DataFrame(universe_list)
        df = df.applymap(str)  # Change from mdanalysis universe type to string
        df[0] = df[0].replace(replace, regex=True)  # replace everything (from '<Residue VAL, 16>' to 'V,16')
        df = pd.concat([df[0].str.split(',', expand=True)], axis=1)  # split column by ',' into two columns
        df = df.drop(df[df[0].map(len) > 1].index)  # remove every none residue row (water, etc)
        df = df.set_index(1)[0].to_dict()  # reforming dataframe and saving as a dict
        df = OrderedDict(df)

        self.pdb_sequence_dict = df

    def assign_generic_numbers_to_pdb(self, sequence_dict, pdb_sequence_dict):
        """
        Assign generic numbers (for class A GPCRs --> Ballosteros Weinstein Numbering)
        :param sequence_dict: Retrieved through get_sequence_dict
        :param pdb_sequence_dict: Retrieved through get_pdb_sequence_dict
        :return: saved dict in instance attribute "pdb_generic_numbers_dict"
        """
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

    def align_with_generic_numbers(self, instance):
        """
        Start all methods in right order to get assigned generic numbers.
        :param instance: Name of instance
        :return: saved dict in instance attribute "pdb_generic_numbers_dict"
        """
        instance.get_protein_name(self.identifier)
        instance.get_sequence_dict(self.protein_name)
        instance.get_pdb_sequence_dict(self.preferred_chain, self.topology_dir)
        instance.assign_generic_numbers_to_pdb(self.sequence_dict, self.pdb_sequence_dict)


class GetDescriptors:
    """
    Get different descriptor using instance attributes from class GetReceptorMapping.
    """
    def __init__(self):
        pass

    def get_dihedral(self):
        pass

    def get_descriptor_x(self):
        pass


# START MY TIMER
start = timer()

s1pr1 = GetReceptorMapping("/home/gocky/python/GPCR-ML_data/pdbs/3V2Y.pdb", "", "3V2Y")
s1pr1.align_with_generic_numbers(s1pr1)
print(s1pr1.pdb_generic_numbers_dict)
print(s1pr1.pdb_sequence_dict)

s1pr2 = GetReceptorMapping("/mdspace/gocky/S1PR2_project/homology_models_e28_docked/s1pr2.pdb", "", "s1pr2_human")
s1pr2.align_with_generic_numbers(s1pr2)
#print(s1pr2.pdb_generic_numbers_dict)


# STOP MY TIMER
elapsed_time = timer() - start  # in seconds
print(elapsed_time)


