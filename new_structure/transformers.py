import exceptions
from abc import ABC, abstractmethod
from itertools import chain, product

from rdkit import Chem

import utils


class IMolTransformer(ABC):
    """
    Interface for classes that transform molecular structure.
    """

    @abstractmethod
    def get_args(self, data):
        """
        Abstract method for formatting the arguments that will be passed to the method transform().

        Args:
            data (iterable): An iterable, possibly containing other iterables, that contain the molecule data as dicts.
        """

    @abstractmethod
    def transform(self, args):
        """
        Abstract method for applying the specific molecular transformation handled by the derived class on the molecules
        passed in as arguments.

        Args:
            args (iterable): An iterable containing the molecules and associated data that are to be operated on.
        """


class SideChainGenerator(IMolTransformer):
    """
    Implementation of an IMolTransformer that takes a connection molecule and heterocycle and creates a set of
    side_chains by attaching the connection molecule to all valid positions on the heterocycle. Valid positions are
    determined in instace method is_valid_atom().
    """

    # atom map numbers
    _CONNECTION_MAP_NUM = 1
    _HETEROCYCLE_MAP_NUM = 2

    def get_args(self, data):

        return product(data.heterocycles, data.connections)

    def transform(self, args):

        heterocycle_doc, connection_doc = args
        heterocycle = Chem.Mol(heterocycle_doc['binary'])
        connection = Chem.Mol(connection_doc['binary'])

        # tag any methyl groups on heterocycle so they are not confused with methyl connection atoms in the monomer
        # generation step
        ignored_map_nums = []
        patt = Chem.MolFromSmarts('[CH3]')
        for i, atom in enumerate(chain.from_iterable(heterocycle.GetSubstructMatches(patt)), start=3):
            ignored_map_nums.append(i)
            heterocycle.GetAtomWithIdx(atom).SetAtomMapNum(i)

        # check if connecting atom is atom mapped
        self.is_valid_connection(connection)

        # make attachment at each atom
        unique_mols = {}
        for atom in heterocycle.GetAtoms():

            # detetmine atom eligibility
            if not self.is_valid_atom(atom):
                continue

            # merge parent side chain with conenction
            atom.SetAtomMapNum(self._HETEROCYCLE_MAP_NUM)
            side_chain = utils.connect_mols(heterocycle, connection, ignored_map_nums=ignored_map_nums)
            atom.SetAtomMapNum(0)

            # check for uniqueness and record results
            binary = side_chain.ToBinary()
            Chem.Kekulize(side_chain)
            unique_mols[binary] = (Chem.MolToSmiles(side_chain, kekuleSmiles=True), atom.GetIdx())

        return self.format_data(unique_mols, heterocycle_doc, connection_doc)

    def is_valid_atom(self, atom):
        """
        Helper method used to determine if a specific position on the heterocycle in valid.

        Args:
            atom (RDKit Atom): The candidate atom on the heterocycle that might have the connection molecule attached
                to it.

        Returns:
            bool: True if valid.
        """

        valid_carbon = atom.GetSymbol() == 'C' and 0 < atom.GetTotalNumHs() < 3
        valid_hetero = atom.GetSymbol() in ('N', 'O', 'S') and atom.GetTotalNumHs() != 0

        return valid_carbon or valid_hetero

    def is_valid_connection(self, connection):
        """
        Helper method that determines whether the connection molecule has an atom map number specifying the point of
        attachment or not.

        Args:
            connection (RDKit Mol): The connection molecule.

        Raises:
            exceptions.MissingMapNumber: Raised when there is no atom map number on the connection molecule.
        """

        if self._CONNECTION_MAP_NUM not in [atom.GetAtomMapNum() for atom in connection.GetAtoms()]:
            raise exceptions.MissingMapNumber('Connection molecule missing atom map number')

    def format_data(self, side_chains, heterocycle, connection):
        """
        Helper method that fills in a dict object for each new side_chain with the necessary data associated with that
        side_chain needed for record keeping.

        Args:
            side_chains (dict): A dict containing all the unquie side_chains and their associated data created from the
                heterocycle and connection molecules.
            heterocycle (dict): The dict containing the data associated with the heterocycle molecule.
            connection (dict): The dict containing the data associated with the connection molecule.

        Returns:
            list: A list containing the newly created side_chain dicts.
        """

        data = []
        # chunk = len(side_chains) * self.connections.index(connection)
        chunk = len(side_chains) * 3
        for i, (binary, (kekule, conn_atom_idx)) in enumerate(side_chains.items()):
            data.append({'_id': heterocycle['_id'] + str(chunk + i),
                         'type': 'side_chain',
                         'binary': binary,
                         'kekule': kekule,
                         'conn_atom_idx': conn_atom_idx,
                         'heterocycle': heterocycle['_id'],
                         'connection': connection['_id']})

        return data


class MonomerGenerator(IMolTransformer):
    def get_args(self, data):
        pass

    def transform(self, args):
        pass


class PeptideGenerator(IMolTransformer):
    def get_args(self, data):
        pass

    def transform(self, args):
        pass


class TemplatePeptideGenerator(IMolTransformer):
    def get_args(self, data):
        pass

    def transform(self, args):
        pass


class MacrocycleGenerator(IMolTransformer):
    def get_args(self, data):
        pass

    def transform(self, args):
        pass


class Methylateor(IMolTransformer):
    def get_args(self, data):
        pass

    def transform(self, args):
        pass


class StereoChemPermutor(IMolTransformer):
    def get_args(self, data):
        pass

    def transform(self, args):
        pass


class MacrocycleConformerGenerator(IMolTransformer):
    def get_args(self, data):
        pass

    def transform(self, args):
        pass
