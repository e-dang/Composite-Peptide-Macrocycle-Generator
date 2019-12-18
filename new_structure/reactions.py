from abc import ABC, abstractmethod
import molecules
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from itertools import chain
import utils
from copy import deepcopy
import exceptions


class IReaction(ABC):

    TEMPLATE_OLIGOMERIZATION_MAP_NUM = molecules.ITemplateMol.OLIGOMERIZATION_MAP_NUM  # 1
    TEMPLATE_EAS_MAP_NUM = molecules.ITemplateMol.EAS_CARBON_MAP_NUM  # 2
    SIDECHAIN_OLIGOMERIZATION_MAP_NUM = 3
    SIDECHAIN_EAS_MAP_NUM = 4

    @abstractmethod
    def __call__(self):
        pass

    @property
    @abstractmethod
    def name(self):
        pass

    @property
    @abstractmethod
    def type(self):
        pass

    @property
    @abstractmethod
    def applicable_template(self):
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

    def __bool__(self):
        return self.valid

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


class AbstractUniMolecularReaction(IReaction):

    CARBON_ALDEHYDE_MAP_NUM = molecules.ITemplateMol.PS_CARBON_ALDEHYDE_MAP_NUM  # 7
    OXYGEN_ALDEHYDE_MAP_NUM = molecules.ITemplateMol.PS_OXYGEN_ALDEHYDE_MAP_NUM  # 8

    @property
    def type(self):
        return 'unimolecular'


class AbstractBiMolecularReaction(IReaction):

    def initialize(self, sidechain, template, reacting_atom):
        self.sidechain = sidechain
        self.template = Chem.MolFromSmiles(template)
        self.reacting_atom = reacting_atom

    @property
    def type(self):
        return 'bimolecular'


class AbstractBMSideChainReaction(AbstractBiMolecularReaction):

    @property
    def type(self):
        return super().type + '_sidechain_reaction'

    def tag_connection_atom(self):

        matches = self.sidechain.GetSubstructMatches(Chem.MolFromSmarts('[CH3]'))
        for atom_idx in chain.from_iterable(matches):
            atom = self.sidechain.GetAtomWithIdx(atom_idx)
            if atom.GetIsotope() == 13:  # isotope labeling for methyl carbons that arent candidate attachments
                atom.SetIsotope(0)
            else:
                atom.SetAtomMapNum(self.SIDECHAIN_OLIGOMERIZATION_MAP_NUM)
                utils.atom_to_wildcard(atom)


class AbstractBMMonomerReaction(AbstractBiMolecularReaction):

    BACKBONE_OLIGOMERIZATION_MAP_NUM = molecules.IBackBoneMol.OLIGOMERIZATION_MAP_NUM  # 1
    BACKBONE_CARBOXYL_MAP_NUM = 5
    BACKBONE_NITROGEN_MAP_NUM = 6

    @abstractmethod
    def create_reactant(self):
        pass

    @property
    def type(self):
        return super().type + '_monomer_reaction'

    def tag_connection_atom(self):

        matches = self.sidechain.GetSubstructMatches(Chem.MolFromSmarts('[CH3]'))
        for atom_idx in chain.from_iterable(matches):
            atom = self.sidechain.GetAtomWithIdx(atom_idx)
            if atom.GetIsotope() == 13:  # isotope labeling for methyl carbons that arent candidate attachments
                atom.SetIsotope(0)
            else:
                atom.SetAtomMapNum(self.SIDECHAIN_OLIGOMERIZATION_MAP_NUM)

    def create_monomer(self):

        # get modified backbone
        backbone = utils.get_partial_backbone(self.BACKBONE_CARBOXYL_MAP_NUM)

        # tag n-terminus nitrogen for oligomerization
        for atom_idx in chain.from_iterable(backbone.GetSubstructMatches(Chem.MolFromSmarts('[NH2]'))):
            atom = backbone.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() == 'N':
                atom.SetAtomMapNum(self.BACKBONE_NITROGEN_MAP_NUM)
                break

        # create monomer
        map_nums = (self.BACKBONE_OLIGOMERIZATION_MAP_NUM, self.SIDECHAIN_OLIGOMERIZATION_MAP_NUM)
        return utils.connect_mols(self.sidechain, backbone, map_nums=map_nums)


class FriedelCrafts(AbstractBMSideChainReaction):

    def __call__(self, sidechain, template, reacting_atom):
        try:
            self.reset()
            super().initialize(sidechain, template.friedel_crafts_kekule, reacting_atom)
            self.validate()
        except AttributeError:  # template doesn't have a friedel_crafts_kekule
            self.valid = False

    @property
    def name(self):
        return 'friedel_crafts'

    @property
    def applicable_template(self):
        return 'all'

    def validate(self):

        # check if sidechain reacting atom is valid
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
        self.tag_connection_atom()
        self.create_product()
        self.create_reaction_smarts([self.sidechain, self.template], self.product)

    def create_product(self):
        map_nums = (self.SIDECHAIN_EAS_MAP_NUM, self.TEMPLATE_EAS_MAP_NUM)
        self.product = utils.connect_mols(self.sidechain, self.template, map_nums=map_nums, clear_map_nums=False)


class TsujiTrost(AbstractBMSideChainReaction):

    def __call__(self, sidechain, template, reacting_atom):
        try:
            self.reset()
            super().initialize(sidechain, template.tsuji_trost_kekule, reacting_atom)
            self.validate()
        except AttributeError:  # doesn't have tsuji_trost_kekule attribute
            self.valid = False

    @property
    def name(self):
        return 'tsuji_trost'

    @property
    def applicable_template(self):
        return 'all'

    def validate(self):

        # check if sidechain reacting atom is valid
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

        self.tag_connection_atom()
        self.create_product()
        self.create_reaction_smarts([self.sidechain, self.template], self.product)

    def create_product(self):

        map_nums = [self.SIDECHAIN_EAS_MAP_NUM, self.TEMPLATE_EAS_MAP_NUM]
        self.product = utils.connect_mols(self.sidechain, self.template, map_nums=map_nums, clear_map_nums=False)


class PictetSpangler(AbstractBMMonomerReaction):

    TEMPLATE_CARBON_ALDEHYDE_MAP_NUM = molecules.ITemplateMol.PS_CARBON_ALDEHYDE_MAP_NUM  # 7
    TEMPLATE_OXYGEN_ALDEHYDE_MAP_NUM = molecules.ITemplateMol.PS_OXYGEN_ALDEHYDE_MAP_NUM  # 8

    def __call__(self, sidechain, template, reacting_atom):
        try:
            self.reset()
            super().initialize(sidechain, template.pictet_spangler_kekule, reacting_atom)
            self.template_name = template.name
            self.tag_connection_atom()
            self.validate()
        except AttributeError:  # doesn't have pictet_spangler_kekule attribute
            self.valid = False

    @property
    def name(self):
        return 'pictet_spangler'

    @property
    def applicable_template(self):
        return self.template_name

    def validate(self):

        # check reacting atom is aromatic carbon with at least one hydrogen
        if not self.reacting_atom.GetIsAromatic() \
                or self.reacting_atom.GetSymbol() != 'C' \
                or self.reacting_atom.GetTotalNumHs() == 0:
            self.valid = False
            return

        # check that attachment point of sidechain is two atoms away from reacting atom
        for atom in self.sidechain.GetAtoms():
            if atom.GetAtomMapNum() == self.SIDECHAIN_OLIGOMERIZATION_MAP_NUM:
                connection_atom = atom
                break
        paths = Chem.FindAllPathsOfLengthN(self.sidechain, 3, useBonds=False, rootedAtAtom=self.reacting_atom.GetIdx())
        atoms = set().union([atom for path in paths for atom in path])
        if connection_atom.GetIdx() not in atoms - set(atom.GetIdx() for atom in connection_atom.GetNeighbors()):
            self.valid = False
            return

        # check that template has unmasked aldehyde (this aldehyde will participate in pictet spangler)
        if not self.template.GetSubstructMatch(Chem.MolFromSmarts('C(=O)')):
            raise exceptions.InvalidMolecule(
                'The provided template must have an aldehyde to participate in a pictet spangler reaction')

        self.valid = True

    def create_reaction(self):
        self.create_reactant()
        self.create_product()
        self.create_reaction_smarts([self.reactant], self.product)

    def create_reactant(self):
        monomer = super().create_monomer()
        map_nums = (self.TEMPLATE_OLIGOMERIZATION_MAP_NUM, self.BACKBONE_NITROGEN_MAP_NUM)
        self.reactant = utils.connect_mols(monomer, self.template, map_nums=map_nums, clear_map_nums=False)

    def create_product(self):

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
        map_nums = (self.SIDECHAIN_EAS_MAP_NUM, self.TEMPLATE_CARBON_ALDEHYDE_MAP_NUM)
        reactant = utils.connect_mols(reactant, map_nums=map_nums, clear_map_nums=False)

        # create bond between peptide nitrogen and unmasked aldehyde carbon
        map_nums = (self.BACKBONE_NITROGEN_MAP_NUM, self.TEMPLATE_CARBON_ALDEHYDE_MAP_NUM)
        self.product = utils.connect_mols(reactant, map_nums=map_nums, clear_map_nums=False)


class TemplatePictetSpangler(AbstractUniMolecularReaction):

    REACTING_ATOM_MAP_NUM = molecules.ITemplateMol.TEMPLATE_PS_REACTING_ATOM_MAP_NUM  # 9
    NITROGEN_MAP_NUM = molecules.ITemplateMol.TEMPLATE_PS_NITROGEN_MAP_NUM  # 10

    def __call__(self, template):
        try:
            self.reset()
            self.template = Chem.MolFromSmiles(template.template_pictet_spangler_kekule)
            self.template_name = template.name
            self.validate()
        except AttributeError:  # doesn't have template_pictet_spangler_kekule attribute
            self.valid = False

    @property
    def name(self):
        return 'template_pictet_spangler'

    @property
    def applicable_template(self):
        return self.template_name

    def validate(self):

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
        self.create_product()
        self.create_reaction_smarts([self.template], self.product)

    def create_product(self):

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

        # create bond between reacting atom and the aldehyde carbon
        map_nums = (self.CARBON_ALDEHYDE_MAP_NUM, self.REACTING_ATOM_MAP_NUM)
        template = utils.connect_mols(template, map_nums=map_nums, clear_map_nums=False)

        # create bond between the peptide nitrogen atom and the aldehyde carbon
        map_nums = (self.NITROGEN_MAP_NUM, self.CARBON_ALDEHYDE_MAP_NUM)
        self.product = utils.connect_mols(template, map_nums=map_nums, clear_map_nums=False)


class UnmaskedAldehydeCyclization(AbstractUniMolecularReaction):

    BACKBONE_CARBOXYL_MAP_NUM = 5
    BACKBONE_NITROGEN_MAP_NUM = 6

    def __call__(self, template):
        try:
            self.reset()
            self.template = Chem.MolFromSmiles(template.pictet_spangler_kekule)
            self.template_name = template.name
            self.validate()
        except AttributeError:  # doesn't have template_pictet_spangler_kekule attribute
            self.valid = False

    @property
    def name(self):
        return 'unmasked_aldehyde_cyclization'

    @property
    def applicable_template(self):
        return self.template_name

    def validate(self):

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
        self.create_reactant()
        self.create_product()
        self.create_reaction_smarts([self.reactant], self.product)

    def create_reactant(self):

        backbone = utils.get_partial_backbone(self.BACKBONE_CARBOXYL_MAP_NUM)

        # tag n-terminus nitrogen for oligomerization and clear backbone carbon oligomerzation map num
        for atom in backbone.GetAtoms():
            if atom.GetSymbol() == 'N':
                atom.SetAtomMapNum(self.BACKBONE_NITROGEN_MAP_NUM)
            elif atom.GetAtomMapNum() == molecules.IBackBoneMol.OLIGOMERIZATION_MAP_NUM:
                atom.SetAtomMapNum(0)

        map_nums = (self.TEMPLATE_OLIGOMERIZATION_MAP_NUM, self.BACKBONE_NITROGEN_MAP_NUM)
        self.reactant = utils.connect_mols(self.template, backbone, map_nums=map_nums, clear_map_nums=False)

    def create_product(self):

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


class PyrroloIndolene(AbstractBMMonomerReaction):

    ADJ_CARBON_MAP_NUM = 7
    N_TERM_WILDCARD_MAP_NUM = 8

    def __call__(self, sidechain, template, reacting_atom):
        try:
            self.reset()
            super().initialize(sidechain, template.pyrrolo_indolene_kekule, reacting_atom)
            self.tag_connection_atom()
            self.validate()
        except AttributeError:  # doesn't have pyrrolo_indolene_kekule attribute
            self.valid = False

    @property
    def name(self):
        return 'pyrrolo_indolene'

    @property
    def applicable_template(self):
        return 'all'

    def validate(self):

        # check template reaction kekule has propylene substructure
        if not self.template.GetSubstructMatch(Chem.MolFromSmarts('C=CC')):
            raise exceptions.InvalidMolecule(
                'The provided template must have a propylene substructure to participate in a pyrrolo-indolene reaction')

        # side_chain must be an indole
        if not self.sidechain.GetSubstructMatch(Chem.MolFromSmarts('*c1c[nH]c2ccccc12')):
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
        for atom in self.sidechain.GetAtoms():
            if atom.GetAtomMapNum() == self.SIDECHAIN_OLIGOMERIZATION_MAP_NUM:
                connection_atom = atom
                break
        if self.reacting_atom.GetIdx() not in [atom.GetIdx() for atom in connection_atom.GetNeighbors()]:
            self.valid = False
            return

        self.valid = True

    def create_reaction(self):

        self.create_reactant()
        self.create_product()
        self.create_reaction_smarts([self.monomer, self.template], self.product)

    def create_reactant(self):

        self.monomer = super().create_monomer()

        # attach an wildcard atom to n-terminus so reaction matches indole at any position in peptide chain
        extra_atom = Chem.MolFromSmarts(f'[CH4:{self.N_TERM_WILDCARD_MAP_NUM}]')  # will be changed to wildcard
        map_nums = (self.BACKBONE_NITROGEN_MAP_NUM, self.N_TERM_WILDCARD_MAP_NUM)
        self.monomer = utils.connect_mols(self.monomer, extra_atom, map_nums=map_nums, clear_map_nums=False)

        for atom in self.monomer.GetAtoms():
            if atom.GetAtomMapNum() == self.N_TERM_WILDCARD_MAP_NUM:  # change methyl to wildcard
                utils.atom_to_wildcard(atom)
            elif atom.GetAtomMapNum() == self.SIDECHAIN_EAS_MAP_NUM:  # tag carbon adjacent to reacting atom
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetTotalNumHs() == 1 and neighbor.GetIsAromatic():
                        neighbor.SetAtomMapNum(self.ADJ_CARBON_MAP_NUM)

    def create_product(self):

        # copy of monomer that will be changed into product
        monomer = deepcopy(self.monomer)

        for atom in monomer.GetAtoms():
            if atom.GetAtomMapNum() == self.SIDECHAIN_EAS_MAP_NUM:
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
        reactant = utils.connect_mols(monomer, map_nums=map_nums, clear_map_nums=False)

        # merge template with monomer
        map_nums = (self.TEMPLATE_EAS_MAP_NUM, self.SIDECHAIN_EAS_MAP_NUM)
        self.product = utils.connect_mols(reactant, self.template, map_nums=map_nums, clear_map_nums=False)
