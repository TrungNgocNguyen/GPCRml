# python standard packages
from collections import OrderedDict
import json
import pickle
from urllib.request import urlopen, urlretrieve

# external packages
from Bio import pairwise2
import MDAnalysis as mda
from MDAnalysis.analysis.dihedrals import Dihedral
import os
import pandas as pd
import seaborn as sns

# GPCRml modules
from .checker import *


class Receptor(object):
    r"""The GPCRml Receptor class contains information describing a GPCR. The Receptor class always requires a path to
    a topology file and the appropriate UniProt ID. The topology file can be of any format supported by MDAnalysis. The
    UniProt ID identifies a receptor and is needed to retrieve necessary information from the GPCRdb API.

    :param topology_path:
        path to topology file
    :type topology_path: ``str``

    :param uniprotid:
        UniProt ID
    :type uniprotid: ``str``

    :param \**kwargs:
        See below

    :Keyword Arguments:
        * *trajectory_path* (``str``) --
            path to trajectory file
        * *pdb_code* (``str``) --
            PDB code
        * *preferred_chain* (``str``) --
            preferred chain for Receptor described by PDB according to GPCRdb
        * *library_directory* (``str``) --
            directory to local receptor library
        * *receptor_name* (``str``) --
            name of receptor according to GPCRdb
        * *reference_sequence* (``str``) --
            amino acid sequence with one-letter codes from the reference structure used for assigning generic numbers
        * *reference_generic_numbers* (``dict``) --
            generic number of reference structure dictionary with residue id as key and generic number as value
        * *topology_sequence* (``str``) --
            amino acid sequence with one-letter codes from the topology structure used for assigning generic numbers
        * *topology_resids* (``list``) --
            residue ids rom the topology structure matching the topology sequence
        * *topology_generic_numbers* (``dict``) --
            generic number of topology structure dictionary with residue id as key and generic number as value
    """
    def __init__(self, topology_path, uniprotid, **kwargs):
        self.topology_path = topology_path
        self.uniprotid = uniprotid
        self.trajectory_path = kwargs.pop('trajectory_path', '')
        self.pdb_code = kwargs.pop('pdb_code', '').upper()
        self.preferred_chain = kwargs.pop('preferred_chain', '')
        self.library_directory = kwargs.pop('library_directory', '')
        self.receptor_name = kwargs.pop('receptor_name', '')
        self.reference_sequence = kwargs.pop('reference_sequence', '')
        self.reference_generic_numbers = kwargs.pop('reference_sequence_dict', None)
        self.topology_sequence = kwargs.pop('topology_sequence', '')
        self.topology_resids = kwargs.pop('topology_resids', None)
        self.topology_generic_numbers = kwargs.pop('topology_sequence_dict', None)

    def __repr__(self):
        return '<Receptor with UniProt ID {}>'.format(self.uniprotid)

    def get_protein_data(self):
        r"""Retrieves receptor name and reference sequence via GPCRdb API according to specified UniProt ID. Generates
        generic numbers dictionary for reference structure.
        """
        response = urlopen('https://gpcrdb.org/services/protein/accession/{}/'.format(self.uniprotid))
        protein = json.loads(response.read().decode('utf-8'))
        self.receptor_name = protein['entry_name']
        self.reference_sequence = protein['sequence']
        response = urlopen('https://gpcrdb.org/services/residues/{}/'.format(self.receptor_name))
        protein = json.loads(response.read().decode('utf-8'))
        self.reference_generic_numbers = {}
        for residue in protein:
            generic_number = residue['display_generic_number']
            if generic_number is not None and 'TM' in residue['protein_segment']:
                generic_number = generic_number.split('.')[0] + 'x' + generic_number.split('x')[1]
            else:
                generic_number = None
            self.reference_generic_numbers[residue['sequence_number']] = generic_number

    def get_topology_data(self):
        r"""Retrieves topology sequence from specified topology_path. If a PDB code is known to the class, receptor
        name, preferred chain, UniProt ID are retrieved via GPCRdb API. The topology file is downloaded from RCSB and
        the path saved to the class.
        """
        if len(self.pdb_code) == 4 and len(self.library_directory) > 0:
            response = urlopen('https://gpcrdb.org/services/structure/{}/'.format(self.pdb_code))
            protein = json.loads(response.read().decode('utf-8'))
            self.receptor_name = protein['protein']
            self.preferred_chain = protein['preferred_chain']
            response = urlopen('https://gpcrdb.org/services/protein/{}/'.format(self.receptor_name))
            protein = json.loads(response.read().decode('utf-8'))
            self.uniprotid = protein['accession']
            if len(self.topology_path) == 0:
                self.topology_path = '{}/{}.pdb'.format(self.library_directory, self.pdb_code)
                if not os.path.exists('{}/{}.pdb'.format(self.library_directory, self.pdb_code)):
                    print('Downloading {} ...'.format(self.pdb_code))
                    urlretrieve('https://files.rcsb.org/download/{}.pdb'.format(self.pdb_code), self.topology_path)

        replace = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E',
                   'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F',
                   'PRO': 'P', 'SEC': 'U', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}
        u = mda.Universe(self.topology_path)
        u = u.select_atoms('protein')
        if len(self.preferred_chain) > 0:
            u = u.select_atoms('segid {}'.format(self.preferred_chain))
        self.topology_sequence = ''.join([replace[resname] for resname in u.residues.resnames])
        self.topology_resids = list(u.residues.resids)

    def assign_generic_numbers(self, mutation_strict=False):
        r"""Assigns generic numbers to topology structure.

        :param mutation_strict:
            if mutations should be handled strict during generic number assignment (default=False)
        :type mutation_strict: ``bool``
        """
        if len(self.topology_sequence) == 0 or self.topology_resids is None:
            self.get_topology_data()
        if self.reference_generic_numbers is None:
            self.get_protein_data()
        alignment = pairwise2.align.globalxs(self.reference_sequence, self.topology_sequence, -10, 0)[0]
        protein_counter = 1
        topology_counter = 0
        self.topology_generic_numbers = {}
        for protein_residue, topology_residue in zip(alignment[0], alignment[1]):
            if self.reference_generic_numbers[protein_counter] is not None:
                if mutation_strict:  # will not assign generic numbers to e.g. stabilizing mutations
                    if protein_residue == topology_residue:
                        self.topology_generic_numbers[self.topology_resids[topology_counter]] = \
                            self.reference_generic_numbers[protein_counter]
                else:  # might also assign generic numbers to e.g. fusion proteins
                    if protein_residue != '-' and topology_residue != '-':
                        self.topology_generic_numbers[self.topology_resids[topology_counter]] = \
                            self.reference_generic_numbers[protein_counter]
            if protein_residue != '-':
                protein_counter += 1
            if topology_residue != '-':
                topology_counter += 1


class PDBReceptor(Receptor):
    r"""The GPCRml PDBReceptor class contains information describing a GPCR. The PDBReceptor class always requires a
    PDB code also known to the GPCRdb API and the path to the library directory to store topology files retrieved from
    RCSB.

    :param pdb_code:
        PDB code
    :type pdb_code: ``str``

    :param library_directory:
        directory to local receptor library
    :type library_directory: ``str``

    :param \**kwargs:
        See below

    :Keyword Arguments:
        * *topology_path* (``str``) --
            path to topology file
        * *uniprotid* (``str``) --
            UniProt ID
        * *trajectory_path* (``str``) --
            path to trajectory file
        * *preferred_chain* (``str``) --
            preferred chain for Receptor described by PDB according to GPCRdb
        * *library_directory* (``str``) --
            directory to local receptor library
        * *receptor_name* (``str``) --
            name of receptor according to GPCRdb
        * *reference_sequence* (``str``) --
            amino acid sequence with one-letter codes from the reference structure used for assigning generic numbers
        * *reference_generic_numbers* (``dict``) --
            generic number of reference structure dictionary with residue id as key and generic number as value
        * *topology_sequence* (``str``) --
            amino acid sequence with one-letter codes from the topology structure used for assigning generic numbers
        * *topology_resids* (``list``) --
            residue ids rom the topology structure matching the topology sequence
        * *topology_generic_numbers* (``dict``) --
            generic number of topology structure dictionary with residue id as key and generic number as value
    """
    def __init__(self, pdb_code, library_directory, topology_path='', uniprotid='', **kwargs):
        Receptor.__init__(self, topology_path='', uniprotid='', **kwargs)
        self.pdb_code = pdb_code.upper()
        self.library_directory = library_directory
        self.topology_path = topology_path
        self.uniprotid = uniprotid
        self.trajectory_path = kwargs.pop('trajectory_path', '')
        self.preferred_chain = kwargs.pop('preferred_chain', '')
        self.receptor_name = kwargs.pop('receptor_name', '')
        self.reference_sequence = kwargs.pop('reference_sequence', '')
        self.reference_sequence_dict = kwargs.pop('reference_sequence_dict', None)
        self.topology_sequence = kwargs.pop('topology_sequence', '')
        self.topology_resids = kwargs.pop('topology_resids', None)
        self.topology_sequence_dict = kwargs.pop('topology_sequence_dict', None)

    def __repr__(self):
        return '<Receptor with PDB code {}>'.format(self.pdb_code)


class ReceptorLibrary(object):
    r"""The GPCRml ReceptorLibrary stores instances of Receptor and PDBReceptor objects in self.receptors. It can be
    used to add, save and load receptors.

    :param \**kwargs:
        See below

    :Keyword Arguments:
        * *receptors* (``list``) --
            list of receptors
    """
    def __init__(self, **kwargs):
        self.receptors = kwargs.pop('receptors', None)

    def __repr__(self):
        if self.receptors is None:
            return '<Library with 0 receptors>'
        else:
            return '<Library with {} receptors>'.format(len(self.receptors))

    @staticmethod
    def check_available_receptors():
        r"""Sends a request to the GPCRdb API for available GPCRs in RCSB and returns the PDB codes as list of strings.

        :returns:
            list of PDB codes
        :rtype: ``list``
        """
        response = urlopen('https://gpcrdb.org/services/structure/')
        return [receptor['pdb_code'] for receptor in json.loads(response.read().decode('utf-8'))]

    def add_receptor(self, topology_path, uniprotid, **kwargs):
        r"""Adds a Receptor instance with assigned generic numbers to self.receptors.

        :param topology_path:
            path to topology file
        :type topology_path: ``str``

        :param uniprotid:
            UniProt ID
        :type uniprotid: ``str``

        :param \**kwargs:
            same as Receptor class
        """
        receptor = Receptor(topology_path, uniprotid, **kwargs)
        receptor.assign_generic_numbers()
        if self.receptors is None:
            self.receptors = [receptor]
        else:
            self.receptors.append(receptor)
        print('Added receptor with UniProt ID {} to receptor library.'.format(uniprotid))

    def add_receptors_from_pdb(self, pdb_codes, library_directory, **kwargs):
        r"""Adds PDBReceptor instances with assigned generic numbers to self.receptors based on specified PDB codes.
        A library directory needs to specified to save PDB files retrieved from RCSB.

        :param pdb_codes:
            list of PDB codes
        :type pdb_codes: ``list``

        :param library_directory:
            directory to local receptor library
        :type library_directory: ``str``

        :param \**kwargs:
            same as Receptor class
        """
        for pdb_code in pdb_codes:
            pdb_receptor = PDBReceptor(pdb_code, library_directory, **kwargs)
            pdb_receptor.assign_generic_numbers()
            if self.receptors is None:
                self.receptors = [pdb_receptor]
            else:
                self.receptors.append(pdb_receptor)
            print('Added receptor with PDB code {} to receptor library.'.format(pdb_code.upper()))

    def save_receptors(self, library_path):
        r"""Saves self.receptors as pickled file whose path is described by the specified library path.

        :param library_path:
            path to local receptor library
        :type library_path: ``str``
        """
        pickle.dump(self.receptors, open(library_path, 'wb'))

    def load_receptors(self, library_path):
        r"""Loads receptors from specified library path. The file should be a pickled file containing a list of
        Receptor and PDBReceptor instances.

        :param library_path:
            path to local receptor library
        :type library_path: ``str``
        """
        self.receptors = pickle.load(open(library_path, 'rb'))


class Descriptors:
    r"""The GPCRml Descriptors class allows for retrieval of descriptors mapped to generic numbers. The Descriptors
    class always requires a ReceptorLibrary instance.

    :param receptor_library:
        ReceptorLibrary instance
    :type receptor_library: ``ReceptorLibrary``
    """
    def __init__(self, receptor_library):
        self.receptor_library = receptor_library

    def get_dihedrals(self, dihedral_type, remove_nan=True):
        """Returns dihedral angles for receptors in the self.receptor_library assigned to generic numbers in a pandas
        DataFrame. This functions requires the specification of dihedral type [phi, psi] and if NaNs should be removed
        from the DataFrame (default=True).

        :param dihedral_type:
            dihedral type [phi or psi]
        :type dihedral_type: ``str``

        :param remove_nan:
            specifies if rows with NaNs should be removed (default=True)
        :type remove_nan: ``bool``

        :returns:
            pandas DataFrame with generic numbers as row labels and columns with dihedrals for every receptor and
            every trajectory frame in the receptor library
        :rtype: ``pd.DataFrame``
        """
        dfs = []
        for receptor in self.receptor_library.receptors:
            name = receptor.topology_path.split('/')[-1].split('.')[0]
            if len(receptor.trajectory_path) == 0:
                u = mda.Universe(receptor.topology_path)
            else:
                u = mda.Universe(receptor.topology_path, receptor.trajectory_path)
            protein = u.select_atoms('protein')
            if len(receptor.preferred_chain) > 0:
                protein = protein.select_atoms('segid {}'.format(receptor.preferred_chain))
            resids = [resid for resid, generic_number in sorted(receptor.topology_generic_numbers.items(),
                                                                key=lambda t: t[1])]
            generic_numbers = []
            atom_groups = []
            for resid in resids:
                if dihedral_type == 'phi':
                    atom_group = protein.select_atoms('resid {}'.format(resid)).residues[0].phi_selection()
                elif dihedral_type == 'psi':
                    atom_group = protein.select_atoms('resid {}'.format(resid)).residues[0].psi_selection()
                else:
                    print('Only phi and psi are valid values for dihedral type.')
                    return
                if atom_group is not None:
                    atom_groups.append(atom_group)
                    generic_numbers.append(receptor.topology_generic_numbers[resid])
            dihedrals = Dihedral(atom_groups).run()
            df = pd.DataFrame(data=dihedrals.angles.transpose(),
                              columns=['{}_{}'.format(name, frame) for frame in range(len(u.trajectory))],
                              index=generic_numbers)
            dfs.append(df)
        merged_df = dfs[0]
        for df in dfs[1:]:
            merged_df = merged_df.join(df, how='outer')
        if remove_nan:
            return merged_df.dropna()
        else:
            return merged_df


class GetReceptorMapping:
    """Map PDB sequence to BWN numbering scheme from GPCRdb"""

    def __init__(self, identifier, topology_dir, trajectory_dir=None):
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
            if residue['display_generic_number'] is not None:
                dictionary[residue['sequence_number']] = [residue['amino_acid'],
                                                          residue['display_generic_number'].split('.')[0] + 'x' +
                                                          residue['display_generic_number'].split('x')[1]]
            else:
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


class GetDescriptors(GetReceptorMapping):
    """
    Get different descriptor using instance attributes from class GetReceptorMapping.
    """
    def __init__(self, getreceptormapping, bwn_phi, bwn_psi):
        self.pdb_generic_numbers_dict = getreceptormapping.pdb_generic_numbers_dict
        self.topology_dir = getreceptormapping.topology_dir
        self.trajectory_dir = getreceptormapping.trajectory_dir
        self.pdb_marker = getreceptormapping.pdb_marker
        self.preferred_chain = getreceptormapping.preferred_chain

        self.bwn_phi = bwn_phi
        self.bwn_psi = bwn_psi

        self.dihedrals_phi = None
        self.dihedrals_psi = None

    def get_dihedrals_phi(self):
        generic_dict = self.pdb_generic_numbers_dict
        key_values_swap = {v: k for k, v in generic_dict.items()}
        positions = [int(key_values_swap[x]) for x in self.bwn_phi]

        if self.trajectory_dir is None:
            u = mda.Universe(self.topology_dir)
        else:
            u = mda.Universe(self.topology_dir, self.trajectory_dir)

        ags = []
        for res in u.residues:
            if self.pdb_marker:
                if res.segid == self.preferred_chain and res.resid in positions:
                    ags.append(res.phi_selection())
            else:
                if res.resid in positions:
                    ags.append(res.phi_selection())

        R = Dihedral(ags).run()
        df = pd.DataFrame(R.angles)
        df.columns = self.bwn_phi
        self.dihedrals_phi = df

    def map_dihedrals_phi(self):
        sns.set(style="whitegrid")

        values = self.dihedrals_phi
        data = pd.DataFrame(values, columns=self.bwn_phi)
        data = data.rolling(1).mean()

        sns.lineplot(data=data, palette="tab10", linewidth=0.5)

    def get_dihedrals_psi(self):
        generic_dict = self.pdb_generic_numbers_dict
        key_values_swap = {v: k for k, v in generic_dict.items()}
        positions = [int(key_values_swap[x]) for x in self.bwn_psi]

        if self.trajectory_dir is None:
            u = mda.Universe(self.topology_dir)
        else:
            u = mda.Universe(self.topology_dir, self.trajectory_dir)

        ags = []
        for res in u.residues:
            if self.pdb_marker:
                if res.segid == self.preferred_chain and res.resid in positions:
                    ags.append(res.psi_selection())
            else:
                if res.resid in positions:
                    ags.append(res.psi_selection())

        R = Dihedral(ags).run()
        df = pd.DataFrame(R.angles)
        df.columns = self.bwn_psi
        self.dihedrals_psi = df

    def map_dihedrals_psi(self):
        sns.set(style="whitegrid")

        values = self.dihedrals_psi
        data = pd.DataFrame(values, columns=self.bwn_psi)
        data = data.rolling(1).mean()

        sns.lineplot(data=data, palette="tab10", linewidth=0.5)
