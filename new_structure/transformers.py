import exceptions
from abc import ABC, abstractmethod
from itertools import chain, product

from rdkit import Chem

import utils
import molecules


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
    _CONNECTION_MAP_NUM = molecules.IConnectionMol.OLIGOMERIZATION_MAP_NUM
    _HETEROCYCLE_MAP_NUM = 2

    def get_args(self, data):

        return product(data.heterocycles, data.connections)

    def transform(self, args):

        self.sidechains = {}
        self.heterocycle, self.connection = args
        heterocycle = Chem.Mol(self.heterocycle['binary'])
        connection = self.connection.mol

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
                         'connection': self.connection.name})

        return data


class MonomerGenerator(IMolTransformer):
    """
    Implementation of an IMolTransformer that takes a sidechain and backbone molecule and creates a monomer
    by attaching the alkyl portion of the sidechain to the designated position on the backbone molecule.
    """

    _BB_MAP_NUM = molecules.IBackBoneMol.OLIGOMERIZATION_MAP_NUM
    _SC_MAP_NUM = 2

    def get_args(self, data):

        return product(data.sidechains, data.backbones)

    def transform(self, args):

        self.sidechain, self.backbone = args
        sidechain = Chem.Mol(self.sidechain['binary'])
        backbone = self.backbone.mol

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
        self.monomer = utils.connect_mols(sidechain, backbone)
        atom.SetAtomMapNum(0)

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
        binary = self.monomer.ToBinary()
        Chem.Kekulize(self.monomer)
        return [{'_id': bb_ind + self.sidechain['_id'].upper(),
                 'type': 'monomer',
                 'binary': binary,
                 'kekule': Chem.MolToSmiles(self.monomer, kekuleSmiles=True),
                 'required': bool(list(self.monomer.GetAromaticAtoms())),
                 'backbone': self.backbone.name,
                 'sidechain': {'_id': self.sidechain['_id'],
                               'heterocycle': self.sidechain['heterocycle'],
                               'conn_atom_idx': self.sidechain['conn_atom_idx']}}]


class PeptideGenerator(IMolTransformer):

    _BACKBONES = molecules.get_hashed_backbones()
    _MONOMER_NITROGEN_MAP_NUM = 1
    _PEPTIDE_CARBON_MAP_NUM = 2

    def get_args(self, data):
        # temporary for testing
        for i, arg in enumerate(filter(self.is_valid_monomers, product(data.monomers, repeat=data.peptide_length))):
            if i == 100:
                break
            else:
                yield arg

    def transform(self, args):

        self.monomers = args

        # begin conneting each monomer in monomers
        for i, monomer_doc in enumerate(self.monomers):

            self.monomer = Chem.Mol(monomer_doc['binary'])
            self.backbone = self._BACKBONES[monomer_doc['backbone']]

            # start peptide with first monomer
            if i == 0:
                self.pep_size = 1
                self.peptide = self.monomer
                self.backbone_prev = self.backbone
                continue

            # assign atom map numbers
            monomer_old_attach = self.tag_monomer_n_term()
            carboxyl_atom, pep_old_attach = self.tag_peptide_c_term()

            # remove oxygen atom from carboxyl
            self.peptide = Chem.RWMol(self.peptide)
            self.peptide.RemoveAtom(carboxyl_atom)

            # connect peptide and monomer
            self.peptide = utils.connect_mols(self.peptide, self.monomer)
            self.pep_size += 1
            pep_old_attach.SetAtomMapNum(0)
            monomer_old_attach.SetAtomMapNum(0)
            self.backbone_prev = self.backbone

        return self.format_data()

    def is_valid_monomers(self, monomers):
        """
        Determines the validity of the peptide.

        Args:
            monomers (list): Contains the monomer documents.

        Returns:
            bool: True if at least one monomer is required.
        """

        for monomer in monomers:
            if monomer['required']:
                return True

        return False

    def tag_monomer_n_term(self):

        for atom_idx in self.monomer.GetSubstructMatch(self.backbone):
            atom = self.monomer.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() == 'N' and atom.GetTotalNumHs() != 0:
                atom.SetAtomMapNum(self._MONOMER_NITROGEN_MAP_NUM)
                return atom

    def tag_peptide_c_term(self):

        flag = False
        for pair in self.peptide.GetSubstructMatches(self.backbone_prev):
            for atom_idx in pair:
                atom = self.peptide.GetAtomWithIdx(atom_idx)
                if atom.GetSymbol() == 'O' and atom.GetTotalNumHs() == 1:  # finds protonated oxygen on carboxyl group
                    carboxyl_atom = atom_idx
                elif atom.GetSymbol() == 'C' and atom.GetTotalNumHs() == 0 and \
                        atom.GetHybridization() == Chem.HybridizationType.SP2:  # finds carboxyl carbon
                    attachment_atom = atom

                # this conditional statements determines if the carboxyl oxygen and carbon found above are valid
                # since the nitrogen on a peptide chain should only have 0 or 1 hydrogens (0 if monomer is a proline)
                # or have 2 hydrogens if it is the first monomer in the peptide
                elif (self.pep_size > 1 and atom.GetSymbol() == 'N' and atom.GetTotalNumHs() <= 1) \
                        or self.pep_size == 1:
                    flag = True

            if flag:
                attachment_atom.SetAtomMapNum(self._PEPTIDE_CARBON_MAP_NUM)
                break

        return carboxyl_atom, attachment_atom

    def format_data(self):
        monomer_data = [{key: value for key, value in monomer.items() if key in ('_id', 'side_chain')}
                        for monomer in self.monomers]
        pep_id = ''.join([monomer['_id'] for monomer in monomer_data])

        binary = self.peptide.ToBinary()
        Chem.Kekulize(self.peptide)
        return [{'_id': pep_id,
                 'type': 'peptide',
                 'binary': binary,
                 'kekule': Chem.MolToSmiles(self.peptide, kekuleSmiles=True),
                 'monomers': monomer_data}]


class TemplatePeptideGenerator(IMolTransformer):

    # any primary amine or proline n-terminus
    _ELIGIBLE_NITROGENS = Chem.MolFromSmarts('[$([NH2]),$([NH;R]);!$([NH2]C(=O)*)]')
    _PEPTIDE_NITROGEN_MAP_NUM = 2

    def get_args(self, data):
        return product(data.peptides, data.templates)

    def transform(self, args):

        self.template_peptides = {}
        self.peptide, self.template = args
        peptide = Chem.Mol(self.peptide['binary'])
        template = self.template.oligomerization_mol

        # for each eligible nitrogen form a connection
        for atom_idx in chain.from_iterable(peptide.GetSubstructMatches(self._ELIGIBLE_NITROGENS)):
            atom = peptide.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() == 'N':  # primary amines only
                atom.SetAtomMapNum(self._PEPTIDE_NITROGEN_MAP_NUM)
                break

        # combine and record results
        template_peptide = utils.connect_mols(template, peptide)
        atom.SetAtomMapNum(0)
        binary = template_peptide.ToBinary()
        Chem.Kekulize(template_peptide)
        self.template_peptides[binary] = Chem.MolToSmiles(template_peptide, kekuleSmiles=True)

        return self.format_data()

    def format_data(self):
        # chunk = len(tp_hybrid) * self.templates.index(template)
        chunk = 3
        data = []
        for i, (binary, kekule) in enumerate(self.template_peptides.items()):
            data.append({'_id': self.peptide['_id'] + str(chunk + i),
                         'type': 'template_peptide',
                         'binary': binary,
                         'kekule': kekule,
                         'peptide': {'_id': self.peptide['_id'],
                                     'monomers': self.peptide['monomers']},
                         'template': self.template.name})
        return data


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
