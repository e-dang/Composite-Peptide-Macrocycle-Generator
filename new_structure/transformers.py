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
    sidechains by attaching the connection molecule to all valid positions on the heterocycle. Valid positions are
    determined in instace method is_valid_atom().
    """

    # atom map numbers
    _CONNECTION_MAP_NUM = 1
    _HETEROCYCLE_MAP_NUM = 2

    def get_args(self, data):

        return product(data.heterocycles, data.connections)

    def transform(self, args):

        self.sidechains = {}
        self.heterocycle, self.connection = args
        heterocycle = Chem.Mol(self.heterocycle['binary'])
        connection = Chem.Mol(self.connection['binary'])

        # check if connecting atom is atom mapped
        self.is_valid_connection(connection)

        # tag any methyl groups on heterocycle so they are not confused with methyl connection atoms in the monomer
        # generation step
        ignored_map_nums = []
        patt = Chem.MolFromSmarts('[CH3]')
        for i, atom in enumerate(chain.from_iterable(heterocycle.GetSubstructMatches(patt)), start=3):
            ignored_map_nums.append(i)
            heterocycle.GetAtomWithIdx(atom).SetAtomMapNum(i)

        # make attachment at each atom
        for atom in heterocycle.GetAtoms():

            # detetmine atom eligibility
            if not self.is_valid_atom(atom):
                continue

            # merge parent side chain with conenction
            atom.SetAtomMapNum(self._HETEROCYCLE_MAP_NUM)
            sidechain = utils.connect_mols(heterocycle, connection, ignored_map_nums=ignored_map_nums)
            atom.SetAtomMapNum(0)

            # check for uniqueness and record results
            binary = sidechain.ToBinary()
            Chem.Kekulize(sidechain)
            self.sidechains[binary] = (Chem.MolToSmiles(sidechain, kekuleSmiles=True), atom.GetIdx())

        return self.format_data()

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

    def format_data(self):
        """
        Helper method that fills in a dict object for each new sidechain with the necessary data associated with that
        sidechain needed for record keeping.

        Returns:
            list: A list containing the newly created sidechain dicts.
        """

        data = []
        # chunk = len(sidechains) * self.connections.index(connection)
        chunk = len(self.sidechains) * 3
        for i, (binary, (kekule, conn_atom_idx)) in enumerate(self.sidechains.items()):
            data.append({'_id': self.heterocycle['_id'] + str(chunk + i),
                         'type': 'sidechain',
                         'binary': binary,
                         'kekule': kekule,
                         'conn_atom_idx': conn_atom_idx,
                         'heterocycle': self.heterocycle['_id'],
                         'connection': self.connection['_id']})

        return data


class MonomerGenerator(IMolTransformer):
    """
    Implementation of an IMolTransformer that takes a sidechain and backbone molecule and creates a monomer
    by attaching the alkyl portion of the sidechain to the designated position on the backbone molecule.
    """

    _BB_MAP_NUM = 1
    _SC_MAP_NUM = 2

    def get_args(self, data):

        return product(data.sidechains, data.backbones)

    def transform(self, args):

        self.monomers = {}
        self.sidechain, self.backbone = args

        sidechain = Chem.Mol(self.sidechain['binary'])
        backbone = Chem.Mol(self.backbone['binary'])

        # find candidate attachment point(s) at terminal end of alkyl chains
        matches = sidechain.GetSubstructMatches(Chem.MolFromSmarts('[CH3]'))
        if not len(matches):
            raise exceptions.InvalidMolecule(
                f'Sidechain {Chem.MolToSmiles(sidechain)} is missing an attachment point')

        # set atom map number for attachment point of side chain
        for atom_idx in chain.from_iterable(matches):
            atom = sidechain.GetAtomWithIdx(atom_idx)
            if atom.GetAtomMapNum() >= 3:  # map nums assigned to methyl carbons that are not attachment points
                atom.SetAtomMapNum(0)
            else:
                atom.SetAtomMapNum(self._SC_MAP_NUM)

        # connect monomer and backbone
        monomer = utils.connect_mols(sidechain, backbone)
        atom.SetAtomMapNum(0)

        # record data
        binary = monomer.ToBinary()
        Chem.Kekulize(monomer)
        self.monomers[binary] = (Chem.MolToSmiles(monomer, kekuleSmiles=True),
                                 bool(list(monomer.GetAromaticAtoms())))

        return self.format_data()

    def format_data(self):
        """
        Helper method that fills in a dict object for each new monomer with the necessary data associated with that
        monomer needed for record keeping.

        Returns:
            list: A list containing the newly created monomer dicts.
        """

        # bb_ind = self.backbones.index(backbone)
        bb_ind = str(3)

        data = []
        for binary, (kekule, required) in self.monomers.items():
            data.append({'_id': bb_ind + self.sidechain['_id'].upper(),
                         'type': 'monomer',
                         'binary': binary,
                         'kekule': kekule,
                         'required': required,
                         'backbone': self.backbone['_id'],
                         'sidechain': {'_id': self.sidechain['_id'],
                                       'heterocycle': self.sidechain['heterocycle'],
                                       'conn_atom_idx': self.sidechain['conn_atom_idx']}})

        return data


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
