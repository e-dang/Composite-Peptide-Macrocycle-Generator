from abc import ABC, abstractmethod
import molecules
from rdkit import Chem
from rdkit.Chem import AllChem
from itertools import chain
import utils
import generators


class IReaction(ABC):

    _SIDECHAIN_EAS_MAP_NUM = generators.JointReactionGenerator.SIDECHAIN_EAS_MAP_NUM

    @abstractmethod
    def __call__(self):
        pass

    @abstractmethod
    def __bool__(self):
        pass

    @property
    @abstractmethod
    def name(self):
        pass

    @abstractmethod
    def validate(self):
        pass

    @abstractmethod
    def create_reaction(self):
        pass

    @abstractmethod
    def create_product(self):
        pass

    @property
    def binary(self):
        return self.reaction.ToBinary()

    def create_reaction_smarts(self, reactants, product):

        reactants = [Chem.MolToSmiles(reactant) for reactant in reactants]
        self.smarts = '(' + '.'.join(reactants) + ')>>' + Chem.MolToSmiles(product)
        self.reaction = AllChem.ReactionFromSmarts(self.smarts)

    def reset(self):
        self.reacting_atom = None
        self.sidechain = None
        self.template = None
        self.product = None
        self.reaction = None
        self.smarts = None
        self.valid = None


class IDisJointReaction(IReaction):

    _TEMPLATE_EAS_MAP_NUM = molecules.ITemplateMol.EAS_CARBON_MAP_NUM

    def initialize(self, sidechain, template, reacting_atom):
        self.template = template
        self.sidechain = sidechain
        self.reacting_atom = reacting_atom

    @property
    def is_joint_reaction(self):
        return False


class FriedelCrafts(IDisJointReaction):

    def __call__(self, sidechain, template, reacting_atom):
        self.reset()
        super().initialize(sidechain, template, reacting_atom)
        self.validate()

    def __bool__(self):
        return self.valid

    @property
    def name(self):
        return 'friedel_crafts'

    def validate(self):

        # check if sidechain reacting atom is valid
        if self.reacting_atom.GetSymbol() != 'C' \
                or self.reacting_atom.GetTotalNumHs() != 0 \
                or self.reacting_atom.GetIsAromatic():
            self.valid = False
            return

        # check template reaction kekule is valid
        if not self.template.GetSubstructMatch(Chem.MolFromSmarts('C=CC')):
            self.valid = False
            return

        self.valid = True

    def create_reaction(self):
        self.create_product()
        self.create_reaction_smarts([self.sidechain, self.template], self.product)

    def create_product(self):
        map_nums = (self._SIDECHAIN_EAS_MAP_NUM, self._TEMPLATE_EAS_MAP_NUM)
        self.product = utils.connect_mols(self.sidechain, self.template, map_nums=map_nums)


class IJointReaction(IReaction):

    _TEMPLATE_OLIGOMERIZATION_MAP_NUM = molecules.ITemplateMol.OLIGOMERIZATION_MAP_NUM
    _SIDECHAIN_OLIGOMERIZATION_NUM = generators.JointReactionGenerator.SIDECHAIN_OLIGOMERIZATION_MAP_NUM
    _BACKBONE_OLIGOMERIZATION_NUM = molecules.IBackBoneMol.OLIGOMERIZATION_MAP_NUM
    _BACKBONE_CARBOXYL_MAP_NUM = 5
    _BACKBONE_NITROGEN_MAP_NUM = 6

    @abstractmethod
    def create_reactant(self):
        pass

    def initialize(self, sidechain, template, reacting_atom):
        self.template = template
        self.sidechain = sidechain
        self.reacting_atom = reacting_atom

    @property
    def is_joint_reaction(self):
        return True

    def create_monomer(self):

        # get and modify backbone
        backbone = molecules.AlphaBackBone().tagged_mol
        carboxyl = Chem.MolFromSmarts('C(=O)O')
        replacement = Chem.MolFromSmarts(f'[*:{self._BACKBONE_CARBOXYL_MAP_NUM}]')
        backbone = AllChem.ReplaceSubstructs(backbone, carboxyl, replacement)[0]

        # tag n-terminus nitrogen for oligomerization
        for atom in chain.from_iterable(backbone.GetSubstructMatches(Chem.MolFromSmarts('[NH2]'))):
            atom = backbone.GetAtomWithIdx(atom)
            if atom.GetSymbol() == 'N':
                atom.SetAtomMapNum(self._BACKBONE_NITROGEN_MAP_NUM)
                break

        # create monomer
        map_nums = (self._BACKBONE_OLIGOMERIZATION_NUM, self._SIDECHAIN_OLIGOMERIZATION_NUM)
        return utils.connect_mols(self.sidechain, backbone, map_nums=map_nums)


class PictetSpangler(IJointReaction):

    _TEMPLATE_CARBON_ALDEHYDE_MAP_NUM = molecules.ITemplateMol.PS_CARBON_ALDEHYDE_MAP_NUM
    _TEMPLATE_OXYGEN_ALDEHYDE_MAP_NUM = molecules.ITemplateMol.PS_OXYGEN_ALDEHYDE_MAP_NUM

    def __call__(self, sidechain, template, reacting_atom):
        self.reset()
        super().initialize(sidechain, template, reacting_atom)
        self.validate()

    def __bool__(self):
        return self.valid

    @property
    def name(self):
        return 'pictet_spangler'

    def validate(self):

        # check reacting atom is aromatic carbon
        if not self.reacting_atom.GetIsAromatic() \
                or self.reacting_atom.GetSymbol() != 'C' \
                or self.reacting_atom.GetTotalNumHs() == 0:
            self.valid = False
            return

        # check that attachment point of sidechain is two atoms away from reacting atom
        for atom in self.sidechain.GetAtoms():
            if atom.GetAtomMapNum() == self._SIDECHAIN_OLIGOMERIZATION_NUM:
                connection_atom = atom
                break
        paths = Chem.FindAllPathsOfLengthN(self.sidechain, 3, useBonds=False, rootedAtAtom=self.reacting_atom.GetIdx())
        atoms = set().union([atom for path in paths for atom in path])
        if connection_atom.GetIdx() not in atoms - set(atom.GetIdx() for atom in connection_atom.GetNeighbors()):
            self.valid = False
            return

        # check that template has unmasked aldehyde (this aldehyde will participate in pictet spangler)
        if not self.template.GetSubstructMatch(Chem.MolFromSmarts('C(=O)')):
            self.valid = False
            return

        self.valid = True

    def create_reaction(self):
        self.create_reactant()
        self.create_product()
        self.create_reaction_smarts([self.reactant], self.product)

    def create_reactant(self):
        monomer = super().create_monomer()
        map_nums = (self._TEMPLATE_OLIGOMERIZATION_MAP_NUM, self._BACKBONE_NITROGEN_MAP_NUM)
        self.reactant = utils.connect_mols(monomer, self.template, map_nums=map_nums, clear_map_nums=False)

    def create_product(self):

        reactant = Chem.RWMol(self.reactant)  # copy of reactant molecule that will be transformed into product

        # reset atom map number on the INSTANCE VARIABLE reactant
        for atom in self.reactant.GetAtoms():
            if atom.GetAtomMapNum() == self._TEMPLATE_OXYGEN_ALDEHYDE_MAP_NUM:
                atom.SetAtomMapNum(0)
                break

        # remove oxygen from unmasked aldehyde on LOCAL VARIABLE reactant
        for atom in list(reactant.GetAtoms()):
            if atom.GetAtomMapNum() == self._TEMPLATE_OXYGEN_ALDEHYDE_MAP_NUM:
                reactant.RemoveAtom(atom.GetIdx())
                break

        # create bond between sidechain EAS carbon and unmasked aldehyde carbon
        map_nums = (self._SIDECHAIN_EAS_MAP_NUM, self._TEMPLATE_CARBON_ALDEHYDE_MAP_NUM)
        reactant = utils.connect_mols(reactant, map_nums=map_nums, clear_map_nums=False)

        # create bond between peptide nitrogen and unmasked aldehyde carbon
        map_nums = (self._BACKBONE_NITROGEN_MAP_NUM, self._TEMPLATE_CARBON_ALDEHYDE_MAP_NUM)
        self.product = utils.connect_mols(reactant, map_nums=map_nums, clear_map_nums=False)
