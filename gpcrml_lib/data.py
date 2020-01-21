# python standard packages
from collections import OrderedDict
import json
import pickle
from urllib.request import urlopen, urlretrieve
import sys

# external packages
from Bio import pairwise2
import MDAnalysis as mda
from MDAnalysis.analysis.dihedrals import Dihedral
import os
import pandas as pd


class Receptor:
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
                    urlretrieve('https://files.rcsb.org/download/{}.pdb'.format(self.pdb_code), self.topology_path)

        replace = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E',
                   'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F',
                   'PRO': 'P', 'SEC': 'U', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}
        u = mda.Universe(self.topology_path)
        u = u.select_atoms('protein')
        if len(self.preferred_chain) > 0:
            u = u.select_atoms('segid {}'.format(self.preferred_chain))

        sequence = []
        for resname in u.residues.resnames:
            try:
                sequence.append(replace[resname])
            except KeyError:
                sequence.append('X')

        self.topology_sequence = ''.join(sequence)
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


class ReceptorLibrary:
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
        for counter, pdb_code in enumerate(pdb_codes):
            text = 'Adding receptor with PDB code {} to receptor library ({}/{}).'.format(pdb_code.upper(), counter + 1,
                                                                                          len(pdb_codes))
            sys.stdout.write('\r' + text)
            sys.stdout.flush()
            pdb_receptor = PDBReceptor(pdb_code, library_directory, **kwargs)
            pdb_receptor.assign_generic_numbers()
            if self.receptors is None:
                self.receptors = [pdb_receptor]
            else:
                self.receptors.append(pdb_receptor)
        sys.stdout.write('\n')

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
        num_receptors = len(self.receptor_library.receptors)
        for counter, receptor in enumerate(self.receptor_library.receptors):
            name = receptor.topology_path.split('/')[-1].split('.')[0]
            text = 'Adding {} angles for {} to receptor library ({}/{}).'.format(dihedral_type, name, counter + 1,
                                                                                 num_receptors)
            sys.stdout.write('\r' + text)
            sys.stdout.flush()
            if len(receptor.trajectory_path) == 0:
                u = mda.Universe(receptor.topology_path)
            else:
                u = mda.Universe(receptor.topology_path, receptor.trajectory_path)
            protein = mda.Merge(u.select_atoms('protein and not altloc B'))  # fix for altloc
            if len(receptor.preferred_chain) > 0:
                protein = protein.select_atoms('segid {}'.format(receptor.preferred_chain))
            resids = [resid for resid, generic_number in sorted(receptor.topology_generic_numbers.items(),
                                                                key=lambda t: t[1])]
            generic_numbers = []
            atom_groups = []
            for resid in resids:
                if dihedral_type == 'phi':
                    if resid > 0:
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
        sys.stdout.write('\nMerging results...\n')
        merged_df = dfs[0]
        for df in dfs[1:]:
            merged_df = merged_df.join(df, how='outer')

        if dihedral_type == 'phi':
            merged_df.index += '_phi'
        else:
            merged_df.index += '_psi'

        if remove_nan:
            return merged_df.dropna().swapaxes("index", "columns")
        else:
            return merged_df.swapaxes("index", "columns")
