from abc import ABC, abstractmethod
from copy import deepcopy
from itertools import chain

from rdkit import Chem
from rdkit.Chem import AllChem

import macrocycles.exceptions as exceptions
import macrocycles.molecules as molecules
import macrocycles.utils as utils


class IReaction(ABC):
    """
    Interface for Reaction classes.
    """

    TEMPLATE_OLIGOMERIZATION_MAP_NUM = molecules.ITemplateMol.OLIGOMERIZATION_MAP_NUM  # 1
    TEMPLATE_EAS_MAP_NUM = molecules.ITemplateMol.EAS_CARBON_MAP_NUM  # 2
    SIDECHAIN_OLIGOMERIZATION_MAP_NUM = 3
    REACTING_MOL_EAS_MAP_NUM = 4
    BACKBONE_CARBOXYL_MAP_NUM = 5
    BACKBONE_NITROGEN_MAP_NUM = 6

    @abstractmethod
    def __call__(self):
        """
        Abstract method to initiate the validation of the reaction.
        """

    @property
    @abstractmethod
    def name(self):
        """
        The name of the reaction.
        """

    @property
    @abstractmethod
    def type(self):
        """
        The classification of the reaction.
        """

    @property
    @abstractmethod
    def applicable_template(self):
        """
        The type of templates that are applicable to the reaction.
        """

    @abstractmethod
    def validate(self):
        """
        Abstract method used for determining whether the given molecules can undergo the reaction.
        """

    @abstractmethod
    def create_reaction(self):
        """
        Abstract method that initiates the creation of the reaction SMARTS string.
        """

    @abstractmethod
    def create_product(self):
        """
        Abstract method that creates the product of the reaction from its reactants.
        """

    def __bool__(self):
        """
        Method for returning the validity of the reaction given the reactants.

        Returns:
            bool: True if valid.
        """

        return self.valid

    @property
    def binary(self):
        """
        Property that returns the RDKit Reaction's binary string.

        Returns:
            binary: The RDKit Reaction as a binary string.
        """

        return self.reaction.ToBinary()

    def create_reaction_smarts(self, reactants, product):
        """
        Method for creating the reaction SMARTS string from the reactants and products.

        Args:
            reactants (iterable[RDKit Mol]): The reactants.
            product (RDKit Mol): The product.
        """

        reactants = [Chem.MolToSmiles(reactant) for reactant in reactants]
        self.smarts = '(' + '.'.join(reactants) + ')>>' + Chem.MolToSmiles(product)
        self.reaction = AllChem.ReactionFromSmarts(self.smarts)

    def reset(self):
        """
        Method for resetting all instance variables back to None.
        """

        self.reacting_atom = None
        self.template = None
        self.product = None
        self.reaction = None
        self.smarts = None
        self.valid = None


class AbstractUniMolecularReaction(IReaction):
    """
    Abstract class for uni-molecular reactions, where uni-molecular is used with the typical chemistry definition, i.e.
    reactions that have only one reactant. However, this will not be the case for bi-molecular reactions...see doc for
    details.
    """

    CARBON_ALDEHYDE_MAP_NUM = molecules.ITemplateMol.PS_CARBON_ALDEHYDE_MAP_NUM  # 7
    OXYGEN_ALDEHYDE_MAP_NUM = molecules.ITemplateMol.PS_OXYGEN_ALDEHYDE_MAP_NUM  # 8

    @property
    def type(self):
        """
        The classification of the reaction.
        """

        return 'unimolecular'


class AbstractBiMolecularReaction(IReaction):
    """
    Abstract class for bi-molecular reactions, where bi-molecular is used to refer to a reaction that takes two
    "reactants" such as a sidechain and the template. However, in reality these reactions are not actually bi-molecular
    since the sidechain and template are part of the same molecule. But for purposes of creating reactions here, we
    denote these reactions as bi-molecular since we must pass two separte RDKit Mols as arguments.
    """

    BACKBONE_OLIGOMERIZATION_MAP_NUM = molecules.IBackBoneMol.OLIGOMERIZATION_MAP_NUM  # 1
    C_TERM_WILDCARD_MAP_NUM = 52
    N_TERM_WILDCARD_MAP_NUM = 53

    def initialize(self, reacting_mol, template, reacting_atom):
        """
        Method for assigning instance variables to each argument.

        Args:
            reacting_mol (RDKit Mol): The reacting molecule that is not the template, i.e. sidechain, monomer, etc...
            template (str): The kekule SMILES string of the template.
            reacting_atom (RDKit Atom): The atom on the reacting molecule that will initiate the reaction with the
                template.
        """

        self.reacting_mol = reacting_mol
        self.template = Chem.MolFromSmiles(template)
        self.reacting_atom = reacting_atom

    @property
    def type(self):
        """
        The classification of the reaction.
        """

        return 'bimolecular'

    def reset(self):
        """
        Method that extends the functionality of the base class' reset() method by also resetting instance variables,
        self.reacting_mol, self.is_monomer, and self.is_sidechain. Then calls base class' reset().
        """

        self.reacting_mol = None
        self.is_monomer = False
        self.is_sidechain = False
        super().reset()

    def determine_reacting_mol_type(self):
        """
        Determines what type of molecule the reacting molecule is, such as sidechain or monomer, and sets the
        corresponding instance variable flag to True.
        """

        backbones = map(Chem.MolFromSmarts, [backbone.kekule for backbone in molecules.get_backbones()])
        for backbone in backbones:
            if self.reacting_mol.HasSubstructMatch(backbone):
                self.is_monomer = True
                return

        self.is_sidechain = True

    def tag_sidechain_connection_atom(self, change_to_wildcard=True):
        """
        Method for tagging the connection atom of a sidechain (which is always a methyl group with a non C13 carbon).
        Optionally, this atom may also be changed to a wildcard if the connection is not to be extended into a monomer.

        Args:
            change_to_wildcard (bool, optional): Whether to change the atom to a wildcard or not. Defaults to True.
        """

        matches = self.reacting_mol.GetSubstructMatches(Chem.MolFromSmarts('[CH3]'))
        for atom_idx in chain.from_iterable(matches):
            atom = self.reacting_mol.GetAtomWithIdx(atom_idx)
            if atom.GetIsotope() == 13:  # isotope labeling for methyl carbons that arent candidate attachments
                atom.SetIsotope(0)
            else:
                atom.SetAtomMapNum(self.SIDECHAIN_OLIGOMERIZATION_MAP_NUM)
                if change_to_wildcard:
                    utils.atom_to_wildcard(atom)

    def tag_monomer_connection_atom(self):
        """
        Method for changing the connection atoms of a monomer to wildcards tagged with atom map numbers defined as class
        constants. The connection points on a monomer are the n- and c- termini.
        """

        n_term_patt1 = Chem.MolFromSmarts('[NH2]')  # regular amino acid backbone type
        n_term_patt2 = Chem.MolFromSmarts('[NH;R]')  # proline backbone type
        c_term_patt = Chem.MolFromSmarts('C(=O)[OH]')
        n_term_replacement1 = Chem.MolFromSmarts(
            f'[NH1:{self.BACKBONE_NITROGEN_MAP_NUM}][*:{self.N_TERM_WILDCARD_MAP_NUM}]')
        n_term_replacement2 = Chem.MolFromSmarts(
            f'[N;R{self.BACKBONE_NITROGEN_MAP_NUM}][*:{self.N_TERM_WILDCARD_MAP_NUM}]')
        c_term_replacement = Chem.MolFromSmarts(
            f'[C:{self.BACKBONE_CARBOXYL_MAP_NUM}](=O)[*:{self.C_TERM_WILDCARD_MAP_NUM}]')

        # if reacting_mol doesnt have the specified pattern, ReplaceSubstracts() returns the mol unmodified
        self.reacting_mol = AllChem.ReplaceSubstructs(self.reacting_mol, n_term_patt1, n_term_replacement1)[0]
        self.reacting_mol = AllChem.ReplaceSubstructs(self.reacting_mol, n_term_patt2, n_term_replacement2)[0]
        self.reacting_mol = AllChem.ReplaceSubstructs(self.reacting_mol, c_term_patt, c_term_replacement)[0]

    def create_monomer(self):
        """
        Method for creating a monomer from a sidechain with the c-terminus changed to a wildcard. This is useful for
        reactions that need information about the positioning of monomer within the peptide chain in order to determine
        whether the reaction is valid or not.

        Returns:
            RDKit Mol: The monomer.
        """

        # get modified backbone
        backbone = molecules.AlphaBackBone().tagged_mol
        carboxyl = Chem.MolFromSmarts('C(=O)O')
        replacement = Chem.MolFromSmarts(f'[*:{self.C_TERM_WILDCARD_MAP_NUM}]')
        backbone = AllChem.ReplaceSubstructs(backbone, carboxyl, replacement)[0]

        # tag n-terminus nitrogen for oligomerization
        for atom_idx in chain.from_iterable(backbone.GetSubstructMatches(Chem.MolFromSmarts('[NH2]'))):
            atom = backbone.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() == 'N':
                atom.SetAtomMapNum(self.BACKBONE_NITROGEN_MAP_NUM)
                break

        # create monomer
        map_nums = (self.BACKBONE_OLIGOMERIZATION_MAP_NUM, self.SIDECHAIN_OLIGOMERIZATION_MAP_NUM)
        return utils.connect_mols(self.reacting_mol, backbone, map_nums=map_nums)


class FriedelCrafts(AbstractBiMolecularReaction):
    """
    An implementation of the AbstractBiMolecularReaction class that generates reactions between an aromatic carbon with
    at least one hydrogen on a sidechain molecule, and the methyl group on the propylene substructure of the template
    molecule. This class does not take into consideration RegioSQM predictions.
    """

    def __call__(self, reacting_mol, template, reacting_atom):
        """
        Method to initialize the reaction with new molecules. Calls reset(), initialize(), and validate().

        Args:
            reacting_mol (RDKit Mol): The reacting molecule, either sidechain or monomer.
            template (ITemplateMol): An implementation of ITemplateMol that has the property "friedel_crafts_kekule".
            reacting_atom (RDKit Atom): The atom on the reacting molecule that will react with the template.
        """

        try:
            self.reset()
            super().initialize(reacting_mol, template.friedel_crafts_kekule, reacting_atom)
            self.validate()
        except AttributeError:  # template doesn't have a friedel_crafts_kekule
            self.valid = False

    @property
    def name(self):
        """
        The name of the reaction.
        """

        return 'friedel_crafts'

    @property
    def applicable_template(self):
        """
        The type of templates that are applicable to the reaction.
        """

        return 'all'

    def validate(self):
        """
        Method to validate that the given reacting molecule and template are compatible for this type of reaction. The
        requirements are that the reacting atom must be a aromatic carbon with at least one hydrogen, and that the
        template molecule must have a propylene substructure.

        Raises:
            exceptions.InvalidMolecule: Raised if the given template doesn't have the propylene substructure.
        """

        # check if reacting atom is valid
        if self.reacting_atom.GetSymbol() != 'C' \
                or self.reacting_atom.GetTotalNumHs() == 0 \
                or not self.reacting_atom.GetIsAromatic():
            self.valid = False
            return

        # check template reaction kekule has propylene substructure
        if not self.template.GetSubstructMatch(Chem.MolFromSmarts('C=CC')):
            raise exceptions.InvalidMolecule(
                'The provided template must have a propylene substructure to participate in a friedel crafts reaction')

        self.valid = True

    def create_reaction(self):
        """
        Method that intitiates the generation of the product from the reactants and the resulting reaction SMARTS
        string. User must ensure that the reaction is valid before calling this method otherwise errors are expected to
        occur.

        Raises:
            exceptions.UnIdentifiedMolType: Raised if the molecule type was not identified in call to
                determine_reacting_mol_type()
        """

        self.determine_reacting_mol_type()
        if self.is_sidechain:
            self.tag_sidechain_connection_atom()
        elif self.is_monomer:
            self.tag_monomer_connection_atom()
        else:
            raise exceptions.UnIdentifiedMolType('Unable to determine whether the reacting molecule '
                                                 f'{Chem.MolToSmiles(self.reacting_mol)} is a sidechain or a monomer')
        self.create_product()
        self.create_reaction_smarts([self.reacting_mol, self.template], self.product)

    def create_product(self):
        """
        Method that creates the product from the reactants and stores the result in instance variable self.product.
        Product generation is determined by atom map numbers placed on the reacting atom of the reacting mol and class
        constants of the template molecule.
        """

        map_nums = (self.REACTING_MOL_EAS_MAP_NUM, self.TEMPLATE_EAS_MAP_NUM)
        self.product = utils.connect_mols(self.reacting_mol, self.template, map_nums=map_nums, clear_map_nums=False)


class TsujiTrost(AbstractBiMolecularReaction):
    """
    An implementation of the AbstractBiMolecularReaction class that generates reactions between a heteroatom with at
    least one hydrogen on a sidechain or monomer and the propylene substructure on a template molecule. This class does
    not take into consideration pKa predictions.
    """

    def __call__(self, reacting_mol, template, reacting_atom):
        """
        Method to initialize the reaction with new molecules. Calls reset(), initialize(), and validate().

        Args:
            reacting_mol (RDKit Mol): The reacting molecule, either sidechain or monomer.
            template (ITemplateMol): An implementation of ITemplateMol that has the property "tsuji_trost_kekule".
            reacting_atom (RDKit Atom): The atom on the reacting molecule that will react with the template.
        """

        try:
            self.reset()
            super().initialize(reacting_mol, template.tsuji_trost_kekule, reacting_atom)
            self.validate()
        except AttributeError:  # doesn't have tsuji_trost_kekule attribute
            self.valid = False

    @property
    def name(self):
        """
        The name of the reaction.
        """

        return 'tsuji_trost'

    @property
    def applicable_template(self):
        """
        The type of templates that are applicable to the reaction.
        """

        return 'all'

    def validate(self):
        """
        Method to validate that the given reacting molecule and template are compatible for this type of reaction. The
        requirements are that the reacting atom must be a heteroatom with at least one hydrogen, and that the
        template molecule must have a propylene substructure.

        Raises:
            exceptions.InvalidMolecule: Raised if the given template doesn't have the propylene substructure.
        """

        # check if reacting atom is valid
        if self.reacting_atom.GetSymbol() not in ['N', 'O', 'S'] \
                or self.reacting_atom.GetTotalNumHs() == 0:
            self.valid = False
            return

        # check template reaction kekule has propylene substructure
        if not self.template.GetSubstructMatch(Chem.MolFromSmarts('C=CC')):
            raise exceptions.InvalidMolecule(
                'The provided template must have a propylene substructure to participate in a tsuji trost reaction')

        self.valid = True

    def create_reaction(self):
        """
        Method that intitiates the generation of the product from the reactants and the resulting reaction SMARTS
        string. User must ensure that the reaction is valid before calling this method otherwise errors are expected to
        occur.

        Raises:
            exceptions.UnIdentifiedMolType: Raised if the molecule type was not identified in call to
                determine_reacting_mol_type()
        """

        self.determine_reacting_mol_type()
        if self.is_sidechain:
            self.tag_sidechain_connection_atom()
        elif self.is_monomer:
            self.tag_monomer_connection_atom()
        else:
            raise exceptions.UnIdentifiedMolType('Unable to determine whether the reacting molecule '
                                                 f'{Chem.MolToSmiles(self.reacting_mol)} is a sidechain or a monomer')
        self.create_product()
        self.create_reaction_smarts([self.reacting_mol, self.template], self.product)

    def create_product(self):
        """
        Method that creates the product from the reactants and stores the result in instance variable self.product.
        Product generation is determined by atom map numbers placed on the reacting atom of the reacting mol and class
        constants of the template molecule.
        """

        map_nums = [self.REACTING_MOL_EAS_MAP_NUM, self.TEMPLATE_EAS_MAP_NUM]
        self.product = utils.connect_mols(self.reacting_mol, self.template, map_nums=map_nums, clear_map_nums=False)


class PictetSpangler(AbstractBiMolecularReaction):
    """
    Implementation of the AbstractBiMolecularReaction class that generates reactions initiated at an aromatic carbon
    with at least one hydrogen that is six atoms away from the unmasked aldehyde on the template molecule. Not all
    templates are compatible with this reaction as they need to have an unmasked aldehyde to work.
    """

    TEMPLATE_CARBON_ALDEHYDE_MAP_NUM = molecules.ITemplateMol.PS_CARBON_ALDEHYDE_MAP_NUM  # 7
    TEMPLATE_OXYGEN_ALDEHYDE_MAP_NUM = molecules.ITemplateMol.PS_OXYGEN_ALDEHYDE_MAP_NUM  # 8

    def __call__(self, reacting_mol, template, reacting_atom):
        """
        Method to initialize the reaction with new molecules. Sets instance variable self.template_name with the
        corresponding template name and calls reset(), initialize(), determine_reacting_mol_type(), and validate().

        Args:
            reacting_mol (RDKit Mol): The reacting molecule, either sidechain or monomer.
            template (ITemplateMol): An implementation of ITemplateMol that has the property "pictet_spangler_kekule".
            reacting_atom (RDKit Atom): The atom on the reacting molecule that will react with the template.
        """

        try:
            self.reset()
            super().initialize(reacting_mol, template.pictet_spangler_kekule, reacting_atom)
            self.template_name = template.name
            self.determine_reacting_mol_type()
            self.validate()
        except AttributeError:  # doesn't have pictet_spangler_kekule attribute
            self.valid = False

    @property
    def name(self):
        """
        The name of the template.
        """

        return 'pictet_spangler'

    @property
    def applicable_template(self):
        """
        The type of templates that are applicable to the reaction.
        """

        return self.template_name

    def validate(self):
        """
        Method to validate that the given reacting molecule and template are compatible for this type of reaction. The
        requirements are that the reacting atom must be an aromatic carbon with at least one hydrogen, and be six atoms
        away from the template molecule's unmasked aldehyde.

        Raises:
            exceptions.InvalidMolecule: Raised if the given template doesn't have the propylene substructure.
        """

        # check reacting atom is aromatic carbon with at least one hydrogen
        if not self.reacting_atom.GetIsAromatic() \
                or self.reacting_atom.GetSymbol() != 'C' \
                or self.reacting_atom.GetTotalNumHs() == 0:
            self.valid = False
            return

        if self.is_sidechain:  # check that attachment point of sidechain is two atoms away from reacting atom
            self.tag_sidechain_connection_atom(change_to_wildcard=False)
            connection_atom = utils.find_atom(self.reacting_mol, self.SIDECHAIN_OLIGOMERIZATION_MAP_NUM)
            paths = Chem.FindAllPathsOfLengthN(self.reacting_mol, 3, useBonds=False,
                                               rootedAtAtom=self.reacting_atom.GetIdx())
            atoms = set().union([atom for path in paths for atom in path])
            if connection_atom.GetIdx() not in atoms - set(atom.GetIdx() for atom in connection_atom.GetNeighbors()):
                self.valid = False
                return

        if self.is_monomer:  # check that reacting atom is 4 atoms away from n-terminus
            self.tag_monomer_connection_atom()
            n_term_atom = utils.find_atom(self.reacting_mol, self.BACKBONE_NITROGEN_MAP_NUM)
            paths = Chem.FindAllPathsOfLengthN(self.reacting_mol, 5, useBonds=False,
                                               rootedAtAtom=self.reacting_atom.GetIdx())
            atoms = set().union([atom for path in paths for atom in path])
            if n_term_atom.GetIdx() not in atoms - set(atom.GetIdx() for atom in n_term_atom.GetNeighbors()):
                self.valid = False
                return

        # check that template has unmasked aldehyde (this aldehyde will participate in pictet spangler)
        if not self.template.GetSubstructMatch(Chem.MolFromSmarts('C(=O)')):
            raise exceptions.InvalidMolecule(
                'The provided template must have an aldehyde to participate in a pictet spangler reaction')

        self.valid = True

    def create_reaction(self):
        """
        Method that intitiates the generation of the product from the reactants and the resulting reaction SMARTS
        string. User must ensure that the reaction is valid before calling this method otherwise errors are expected to
        occur.
        """

        self.create_reactant()
        self.create_product()
        self.create_reaction_smarts([self.reactant], self.product)

    def create_reactant(self):
        """
        Method that creates the reactant from the sidechain and template molecules since the reaction requires that
        these two molecules be in certain positions with respect to each other. If the reacting molecule is a sidechain,
        it is first converted to a monomer, which is then connected to the template through amide linkage. If the
        reacting molecule is already a monomer then connection to the template is immediately performed.
        """

        map_nums = (self.TEMPLATE_OLIGOMERIZATION_MAP_NUM, self.BACKBONE_NITROGEN_MAP_NUM)

        if self.is_sidechain:
            monomer = super().create_monomer()
            self.reactant = utils.connect_mols(monomer, self.template, map_nums=map_nums, clear_map_nums=False)
        else:
            # remove wildcard atom attached to n-terminus added in call to tag_monomer_connection_atom()
            wildcard_atom = utils.find_atom(self.reacting_atom, self.N_TERM_WILDCARD_MAP_NUM)
            self.reacting_mol = Chem.RWMol(self.reacting_mol)
            self.reacting_mol.RemoveAtom(wildcard_atom.GetIdx())

            # connect template to monomer
            self.reactant = utils.connect_mols(self.reacting_mol, self.template,
                                               map_nums=map_nums, clear_map_nums=False)

    def create_product(self):
        """
        Method that creates the product from the reactant and stores the result in instance variable self.product.
        Product generation is determined by atom map numbers placed on the reacting atom of the reacting mol and class
        constants of the template molecule.
        """

        reactant = Chem.RWMol(self.reactant)  # copy of reactant molecule that will be transformed into product

        # reset atom map number on the INSTANCE VARIABLE reactant
        for atom in self.reactant.GetAtoms():
            if atom.GetAtomMapNum() == self.TEMPLATE_OXYGEN_ALDEHYDE_MAP_NUM:
                atom.SetAtomMapNum(0)
                break

        # remove oxygen from unmasked aldehyde on LOCAL VARIABLE reactant
        for atom in list(reactant.GetAtoms()):
            if atom.GetAtomMapNum() == self.TEMPLATE_OXYGEN_ALDEHYDE_MAP_NUM:
                reactant.RemoveAtom(atom.GetIdx())
                break

        # create bond between sidechain EAS carbon and unmasked aldehyde carbon
        map_nums = (self.REACTING_MOL_EAS_MAP_NUM, self.TEMPLATE_CARBON_ALDEHYDE_MAP_NUM)
        reactant = utils.connect_mols(reactant, map_nums=map_nums, clear_map_nums=False)

        # create bond between peptide nitrogen and unmasked aldehyde carbon
        map_nums = (self.BACKBONE_NITROGEN_MAP_NUM, self.TEMPLATE_CARBON_ALDEHYDE_MAP_NUM)
        self.product = utils.connect_mols(reactant, map_nums=map_nums, clear_map_nums=False)


class PyrroloIndolene(AbstractBiMolecularReaction):
    """
    Implementation of the AbstractBiMolecularReaction class that generates a reaction between the aromatic carbon atom
    on an indole that has no hydrogens and the template's propylene substructure. This reaction destroys the aromaticity
    of the indole.
    """

    ADJ_CARBON_MAP_NUM = 7

    def __call__(self, sidechain, template, reacting_atom):
        """
        Method to initialize the reaction with new molecules. Calls reset(), initialize(), tag_connection_atom(), and
        validate().

        Args:
            sidechain (RDKit Mol): The candidate sidechain that would react with the template.
            template (ITemplateMol): An implementation of ITemplateMol that has the property "pyrrolo_indolene_kekule".
            reacting_atom (RDKit Atom): The atom on the reacting molecule that will react with the template.
        """

        try:
            self.reset()
            super().initialize(sidechain, template.pyrrolo_indolene_kekule, reacting_atom)
            self.tag_sidechain_connection_atom(change_to_wildcard=False)
            self.validate()
        except AttributeError:  # doesn't have pyrrolo_indolene_kekule attribute
            self.valid = False

    @property
    def name(self):
        """
        The name of the reaction.
        """

        return 'pyrrolo_indolene'

    @property
    def applicable_template(self):
        """
        The type of templates that are applicable to the reaction.
        """

        return 'all'

    def validate(self):
        """
        Method to validate that the given reacting molecule and template are compatible for this type of reaction. The
        requirements are that the sidechain must be like the sidechain of tryptophan and the reacting atom is the carbon
        where that is attached to the methylene joining the indole to the amino acid backbone, and that the template
        molecule must have a propylene substructure.

        Raises:
            exceptions.InvalidMolecule: Raised if the given template doesn't have the propylene substructure.
        """

        # check template reaction kekule has propylene substructure
        if not self.template.GetSubstructMatch(Chem.MolFromSmarts('C=CC')):
            raise exceptions.InvalidMolecule('The provided template must have a propylene substructure to participate '
                                             'in a pyrrolo-indolene reaction')

        # sidechain must be an indole
        sidechain = Chem.MolFromSmiles(Chem.MolToSmiles(self.reacting_mol))
        if not sidechain.GetSubstructMatch(Chem.MolFromSmarts('*c1c[nH]c2ccccc12')):
            self.valid = False
            return

        # reaction initiated at carbon that is the attachment point to amino acid backbone, which contains no hydrogens
        if self.reacting_atom.GetTotalNumHs() != 0 or self.reacting_atom.GetSymbol() != 'C':
            self.valid = False
            return

        # reaction also involves carbon adjacent to reacting_atom which must have a hydrogen
        if 1 not in [atom.GetTotalNumHs() for atom in self.reacting_atom.GetNeighbors()]:
            self.valid = False
            return

        # check if carbon is the correct bridged carbon (adjacent to nitrogen)
        for atom in self.reacting_mol.GetAtoms():
            if atom.GetAtomMapNum() == self.SIDECHAIN_OLIGOMERIZATION_MAP_NUM:
                connection_atom = atom
                break
        if self.reacting_atom.GetIdx() not in [atom.GetIdx() for atom in connection_atom.GetNeighbors()]:
            self.valid = False
            return

        self.valid = True

    def create_reaction(self):
        """
        Method that intitiates the generation of the product from the reactants and the resulting reaction SMARTS
        string. User must ensure that the reaction is valid before calling this method otherwise errors are expected to
        occur.
        """

        self.create_reactant()
        self.create_product()
        self.create_reaction_smarts([self.monomer, self.template], self.product)

    def create_reactant(self):
        """
        Method that creates the reactant from the sidechain and template molecules since the reaction requires that
        these two molecules be in certain positions with respect to each other. The sidechain is converted to a monomer,
        which is then connected to the template through amide linkage.
        """

        self.monomer = super().create_monomer()

        # attach a wildcard atom to n-terminus so reaction matches indole at any position in peptide chain
        extra_atom = Chem.MolFromSmarts(f'[CH4:{self.N_TERM_WILDCARD_MAP_NUM}]')  # will be changed to wildcard
        map_nums = (self.BACKBONE_NITROGEN_MAP_NUM, self.N_TERM_WILDCARD_MAP_NUM)
        self.monomer = utils.connect_mols(self.monomer, extra_atom, map_nums=map_nums, clear_map_nums=False)

        for atom in self.monomer.GetAtoms():
            if atom.GetAtomMapNum() == self.N_TERM_WILDCARD_MAP_NUM:  # change methyl to wildcard
                utils.atom_to_wildcard(atom)
            elif atom.GetAtomMapNum() == self.REACTING_MOL_EAS_MAP_NUM:  # tag carbon adjacent to reacting atom
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetTotalNumHs() == 1 and neighbor.GetIsAromatic():
                        neighbor.SetAtomMapNum(self.ADJ_CARBON_MAP_NUM)

    def create_product(self):
        """
        Method that creates the product from the reactant and stores the result in instance variable self.product.
        Product generation is determined by atom map numbers placed on the reacting atom of the reacting mol and class
        constants of the template molecule.
        """

        # copy of monomer that will be changed into product
        monomer = deepcopy(self.monomer)

        for atom in monomer.GetAtoms():
            if atom.GetAtomMapNum() == self.REACTING_MOL_EAS_MAP_NUM:
                eas_atom = atom
            elif atom.GetAtomMapNum() == self.ADJ_CARBON_MAP_NUM:
                adj_atom = atom

        # bond between eas_atom (i.e. the reacting atom) and adj_atom becomes single bond
        bond = monomer.GetBondBetweenAtoms(eas_atom.GetIdx(), adj_atom.GetIdx())
        bond.SetBondType(Chem.BondType.SINGLE)
        eas_atom.SetNumExplicitHs(eas_atom.GetTotalNumHs() + 1)
        adj_atom.SetNumExplicitHs(adj_atom.GetTotalNumHs() + 1)

        # merge backbone nitrogen to adjacent carbon
        map_nums = (self.BACKBONE_NITROGEN_MAP_NUM, self.ADJ_CARBON_MAP_NUM)
        reactant = utils.connect_mols(monomer, map_nums=map_nums, clear_map_nums=False, stereo='CCW')

        # merge template with monomer
        map_nums = (self.TEMPLATE_EAS_MAP_NUM, self.REACTING_MOL_EAS_MAP_NUM)
        self.product = utils.connect_mols(reactant, self.template, map_nums=map_nums, clear_map_nums=False, stereo='CW')


class TemplatePictetSpangler(AbstractUniMolecularReaction):
    """
    Implementation of the AbstractUniMolecularReaction class that generates a pictet spangler reaction with a single
    template and no sidechain/monomer.
    """

    REACTING_ATOM_MAP_NUM = molecules.ITemplateMol.TEMPLATE_PS_REACTING_ATOM_MAP_NUM  # 9
    NITROGEN_MAP_NUM = molecules.ITemplateMol.TEMPLATE_PS_NITROGEN_MAP_NUM  # 10

    def __call__(self, template):
        """
        Method to initialize the reaction with a new template molecule. Calls reset(), and validate(). The template
        molecule must have a property called "template_pictet_spangler_kekule".

        Args:
            template (ITemplateMol): An implementation of ITemplateMol that has the property
                "template_pictet_spangler_kekule".
        """

        try:
            self.reset()
            self.template = Chem.MolFromSmiles(template.template_pictet_spangler_kekule)
            self.template_name = template.name
            self.validate()
        except AttributeError:  # doesn't have template_pictet_spangler_kekule attribute
            self.valid = False

    @property
    def name(self):
        """
        The name of the reaction.
        """

        return 'template_pictet_spangler'

    @property
    def applicable_template(self):
        """
        The type of templates that are applicable to the reaction.
        """

        return self.template_name

    def validate(self):
        """
        Method to validate that the given template molecule can undergo this type of reaction. The requirements are that
        the template molecule has a reacting aromatic carbon atom, an aldehyde, and an amide nitrogen that are all
        correctly atom mapped.

        Raises:
            exceptions.InvalidMolecule: Raised if the given template doesn't have the properly mapped atoms.
        """

        # check that template mol has all atom map numbers required to make reaction
        c_aldehyde, o_aldehyde, reacting, nitrogen = False, False, False, False
        for atom in self.template.GetAtoms():
            if atom.GetSymbol() == 'C' and atom.GetAtomMapNum() == self.CARBON_ALDEHYDE_MAP_NUM:
                c_aldehyde = True
            elif atom.GetSymbol() == 'O' and atom.GetAtomMapNum() == self.OXYGEN_ALDEHYDE_MAP_NUM:
                o_aldehyde = True
            elif atom.GetAtomMapNum() == self.REACTING_ATOM_MAP_NUM:
                reacting = True
            elif atom.GetSymbol() == 'N' and atom.GetAtomMapNum() == self.NITROGEN_MAP_NUM:
                nitrogen = True

        # if error raised, it means the hard coded template_pictet_spangler_kekule is invalid and should be looked at
        if not c_aldehyde:
            raise exceptions.InvalidMolecule(
                'The provided template molecule needs to have a properly atom mapped aldehyde carbon')
        if not o_aldehyde:
            raise exceptions.InvalidMolecule(
                'The provided template molecule needs to have a properly atom mapped aldehyde oxygen')
        if not reacting:
            raise exceptions.InvalidMolecule(
                'The provided template molecule needs to have a properly atom mapped reacting atom')
        if not nitrogen:
            raise exceptions.InvalidMolecule(
                'The provided template molecule needs to have a properly atom mapped peptide nitrogen atom')

        self.valid = True

    def create_reaction(self):
        """
        Method that intitiates the generation of the product from the reactant and the resulting reaction SMARTS
        string. User must ensure that the reaction is valid before calling this method otherwise errors are expected to
        occur.
        """

        self.create_product()
        self.create_reaction_smarts([self.template], self.product)

    def create_product(self):
        """
        Method that creates the product from the reactant and stores the result in instance variable self.product.
        Product generation is determined by atom map numbers placed on the reacting atoms of template molecule.
        """

        # copy of template that turns into product
        template = Chem.RWMol(self.template)

        # reset atom map number on the INSTANCE VARIABLE template
        for atom in self.template.GetAtoms():
            if atom.GetAtomMapNum() == self.OXYGEN_ALDEHYDE_MAP_NUM:
                atom.SetAtomMapNum(0)
                break

        # remove oxygen from unmasked aldehyde on LOCAL VARIABLE template
        for atom in list(template.GetAtoms()):
            if atom.GetAtomMapNum() == self.OXYGEN_ALDEHYDE_MAP_NUM:
                template.RemoveAtom(atom.GetIdx())
            elif atom.GetAtomMapNum() == self.CARBON_ALDEHYDE_MAP_NUM:
                atom.SetNumExplicitHs(atom.GetTotalNumHs() + 2)

        # create bond between the peptide nitrogen atom and the aldehyde carbon
        map_nums = (self.NITROGEN_MAP_NUM, self.CARBON_ALDEHYDE_MAP_NUM)
        template = utils.connect_mols(template, map_nums=map_nums, clear_map_nums=False)

        # create bond between reacting atom and the aldehyde carbon
        map_nums = (self.CARBON_ALDEHYDE_MAP_NUM, self.REACTING_ATOM_MAP_NUM)
        self.product = utils.connect_mols(template, map_nums=map_nums, stereo='CCW', clear_map_nums=False)


class UnmaskedAldehydeCyclization(AbstractUniMolecularReaction):
    """
    Implementation of the AbstractUniMolecularReaction class that generates a cyclization reaction of the unmasked
    aldehyde on a template molecule.
    """

    def __call__(self, template):
        """
        Method to initialize the reaction with a new template molecule. Calls reset(), and validate(). The template
        molecule must have a property called "unmasked_aldehyde_cyclization_kekule".

        Args:
            template (ITemplateMol): An implementation of ITemplateMol that has the property
                "unmasked_aldehyde_cyclization_kekule".
        """

        try:
            self.reset()
            self.template = Chem.MolFromSmiles(template.unmasked_aldehyde_cyclization_kekule)
            self.template_name = template.name
            self.validate()
        except AttributeError:  # doesn't have unmasked_aldehyde_cyclization_kekule attribute
            self.valid = False

    @property
    def name(self):
        """
        The name of the reaction.
        """

        return 'unmasked_aldehyde_cyclization'

    @property
    def applicable_template(self):
        """
        The type of templates that are applicable to the reaction.
        """

        return self.template_name

    def validate(self):
        """
        Method to validate that the given template molecule can undergo this type of reaction. The requirements are that
        the template molecule has an aldehyde that is correctly atom mapped.

        Raises:
            exceptions.InvalidMolecule: Raised if the given template doesn't have the properly mapped atoms.
        """

        olig_carbon_flag, ps_carbon_flag, ps_oxygen_flag = False, False, False
        aldehyde = Chem.MolFromSmarts('[CH1]=O')
        for atom_idx in chain.from_iterable(self.template.GetSubstructMatches(aldehyde)):
            atom = self.template.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() == 'C' and atom.GetAtomMapNum() == self.TEMPLATE_OLIGOMERIZATION_MAP_NUM:
                olig_carbon_flag = True
            elif atom.GetSymbol() == 'C' and atom.GetAtomMapNum() == self.CARBON_ALDEHYDE_MAP_NUM:
                for neighbor in atom.GetNeighbors():  # must have adjacent methylene
                    if neighbor.GetSymbol() == 'C' and neighbor.GetTotalNumHs() == 2:
                        ps_carbon_flag = True
            elif atom.GetSymbol() == 'O' and atom.GetAtomMapNum() == self.OXYGEN_ALDEHYDE_MAP_NUM:
                ps_oxygen_flag = True

        if not olig_carbon_flag:
            raise exceptions.InvalidMolecule(
                'The provided template molecule needs to have a properly atom mapped aldehyde oligomerization carbon')

        if not ps_carbon_flag:
            self.valid = False
            return

        if not ps_oxygen_flag:
            raise exceptions.InvalidMolecule(
                'The provided template molecule needs to have a properly atom mapped aldehyde oxygen')

        self.valid = True

    def create_reaction(self):
        """
        Method that intitiates the generation of the product from the reactant and the resulting reaction SMARTS
        string. User must ensure that the reaction is valid before calling this method otherwise errors are expected to
        occur.
        """

        self.create_reactant()
        self.create_product()
        self.create_reaction_smarts([self.reactant], self.product)

    def create_reactant(self):
        """
        Method that creates the reactant from the template molecule by adding an amide nitrogen that will be used in
        the reaction.
        """

        backbone = molecules.AlphaBackBone().tagged_mol
        carboxyl = Chem.MolFromSmarts('CC(=O)O')
        replacement = Chem.MolFromSmarts(f'[*:{self.BACKBONE_CARBOXYL_MAP_NUM}]')
        backbone = AllChem.ReplaceSubstructs(backbone, carboxyl, replacement)[0]

        # tag n-terminus nitrogen for oligomerization and clear backbone carbon oligomerzation map num
        for atom in backbone.GetAtoms():
            if atom.GetSymbol() == 'N':
                atom.SetAtomMapNum(self.BACKBONE_NITROGEN_MAP_NUM)
            elif atom.GetAtomMapNum() == molecules.IBackBoneMol.OLIGOMERIZATION_MAP_NUM:
                atom.SetAtomMapNum(0)

        map_nums = (self.TEMPLATE_OLIGOMERIZATION_MAP_NUM, self.BACKBONE_NITROGEN_MAP_NUM)
        self.reactant = utils.connect_mols(self.template, backbone, map_nums=map_nums, clear_map_nums=False)

    def create_product(self):
        """
        Method that creates the product from the reactant and stores the result in instance variable self.product.
        Product generation is determined by atom map numbers placed on the reacting atoms of template molecule.
        """

        # create local copy of reactant to transform into product
        reactant = Chem.RWMol(self.reactant)

        # reset atom map number on aldehyde oxygen on INSTANCE VARIABLE reactant
        for atom in self.reactant.GetAtoms():
            if atom.GetAtomMapNum() == self.OXYGEN_ALDEHYDE_MAP_NUM:
                atom.SetAtomMapNum(0)

        # remove pictet spangler aldehyde oxygen and set the pictet spangler aldehyde carbon's single bond with its
        # neighboring carbon to a double bond in LOCAL copy of reactant
        for atom in list(reactant.GetAtoms()):
            if atom.GetAtomMapNum() == self.OXYGEN_ALDEHYDE_MAP_NUM:
                reactant.RemoveAtom(atom.GetIdx())
            elif atom.GetAtomMapNum() == self.CARBON_ALDEHYDE_MAP_NUM:
                carb_atom = atom
                for bond in atom.GetBonds():
                    if bond.GetBondType() == Chem.BondType.SINGLE:
                        bond.SetBondType(Chem.BondType.DOUBLE)

        # adjust hydrogens
        carb_atom.SetNumExplicitHs(carb_atom.GetTotalNumHs() + 1)

        # connect backbone nitrogen to aldehyde carbon
        map_nums = (self.BACKBONE_NITROGEN_MAP_NUM, self.CARBON_ALDEHYDE_MAP_NUM)
        self.product = utils.connect_mols(reactant, map_nums=map_nums, clear_map_nums=False)


def get_reactions():
    """
    Create and return a list of all reactions.

    Returns:
        list: The reactions.
    """

    return [FriedelCrafts(), TsujiTrost(), PictetSpangler(), TemplatePictetSpangler(), PyrroloIndolene(),
            UnmaskedAldehydeCyclization()]


def get_reactions_of_type(rxn_type):
    """
    Closure for creating a function that returns all reactions of a certain type.

    Args:
        rxn_type (str): The desired reaction type.

    Returns:
        func: A function that when called returns an iterable of the desired reactions.
    """

    def reaction_getter():
        return filter(lambda x: rxn_type in x.type, get_reactions())

    return reaction_getter


get_bimolecular_reactions = get_reactions_of_type('bimolecular')
get_unimolecular_reactions = get_reactions_of_type('unimolecular')
