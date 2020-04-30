from copy import deepcopy
from itertools import chain

from rdkit import Chem
from rdkit.Chem import AllChem

import cpmg.models as models
import macrocycles.utils as utils
from cpmg.exceptions import InvalidMolecule


def create_reaction_smarts(reactants, products):
    smarts = []
    reactants = [Chem.MolToSmiles(reactant) for reactant in reactants]
    for product in products:
        smarts.append('(' + '.'.join(reactants) + ')>>' + Chem.MolToSmiles(product))

    return smarts


class InterMolecularReaction:
    NUCLEOPHILE_EAS_MAP_NUM = 50
    NUCLEOPHILE_WC_MAP_NUM = 51
    BACKBONE_NITROGEN_MAP_NUM = 100
    BACKBONE_CARBOXYL_MAP_NUM = 101
    N_TERM_WILDCARD_MAP_NUM = 102
    C_TERM_WILDCARD_MAP_NUM = 103

    def generate(self, reacting_mol, template, reacting_atom, model):
        self._initialize(reacting_mol, template, model)
        self._validate_reacting_atom(reacting_atom)
        self._validate_template()
        self._create_reactants()
        self._create_product()
        return create_reaction_smarts(self.reactants, self.products)

    def _initialize(self, reacting_mol, template, model):
        self.model = model
        self.reacting_mol = reacting_mol
        self._get_template_mol(template)
        self.reactants = []
        self.products = []

    def _get_template_mol(self, template):
        pass

    def _validate_reacting_atom(self, atom):
        pass

    def _validate_template(self):
        pass

    def _create_reactants(self):
        pass

    def _create_product(self):
        pass

    def _tag_sidechain_connection_atom(self, change_to_wildcard=True):

        matches = self.reacting_mol.GetSubstructMatches(models.METHYL)
        for atom_idx in chain.from_iterable(matches):
            atom = self.reacting_mol.GetAtomWithIdx(atom_idx)
            if atom.GetIsotope() == 13:  # isotope labeling for methyl carbons that arent candidate attachments
                atom.SetIsotope(0)
            else:
                atom.SetAtomMapNum(self.NUCLEOPHILE_WC_MAP_NUM)  # pylint: disable=maybe-no-member
                if change_to_wildcard:
                    utils.atom_to_wildcard(atom)

    def _tag_monomer_connection_atom(self):

        n_term_replacement1 = Chem.MolFromSmarts(
            f'[NH1:{self.BACKBONE_NITROGEN_MAP_NUM}][*:{self.N_TERM_WILDCARD_MAP_NUM}]')
        n_term_replacement2 = Chem.MolFromSmarts(
            f'[N;R{self.BACKBONE_NITROGEN_MAP_NUM}][*:{self.N_TERM_WILDCARD_MAP_NUM}]')
        c_term_replacement = Chem.MolFromSmarts(
            f'[C:{self.BACKBONE_CARBOXYL_MAP_NUM}](=O)[*:{self.C_TERM_WILDCARD_MAP_NUM}]')

        # if nucleophile doesnt have the specified pattern, ReplaceSubstracts() returns the mol unmodified
        self.reacting_mol = AllChem.ReplaceSubstructs(self.reacting_mol, models.N_TERM, n_term_replacement1)[0]
        self.reacting_mol = AllChem.ReplaceSubstructs(self.reacting_mol, models.PROLINE_N_TERM, n_term_replacement2)[0]
        self.reacting_mol = AllChem.ReplaceSubstructs(self.reacting_mol, models.CARBOXYL, c_term_replacement)[0]

    def _create_monomer(self):
        """
        Method for creating a monomer from a sidechain with the c-terminus changed to a wildcard. This is useful for
        reactions that need information about the positioning of monomer within the peptide chain in order to determine
        whether the reaction is valid or not.

        Returns:
            RDKit Mol: The monomer.
        """

        # get modified backbone
        backbone = Chem.MolFromSmarts('N[CH2:1]C(=O)[OH]')
        replacement = Chem.MolFromSmarts(f'[*:{self.C_TERM_WILDCARD_MAP_NUM}]')
        backbone = AllChem.ReplaceSubstructs(backbone, models.CARBOXYL, replacement)[0]
        Chem.SanitizeMol(backbone)

        # tag n-terminus nitrogen for oligomerization
        for atom_idx in chain.from_iterable(backbone.GetSubstructMatches(models.N_TERM)):
            atom = backbone.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() == 'N':
                atom.SetAtomMapNum(self.BACKBONE_NITROGEN_MAP_NUM)
                break

        # create monomer
        map_nums = (self.BACKBONE_OLIGOMERIZATION_MAP_NUM,  # pylint: disable=maybe-no-member
                    self.SIDECHAIN_OLIGOMERIZATION_MAP_NUM)  # pylint: disable=maybe-no-member
        print(Chem.MolToSmiles(self.reacting_mol))
        return utils.connect_mols(self.reacting_mol, backbone, map_nums=map_nums)


class FriedelCrafts(InterMolecularReaction):
    REQUIRED_SUBSTRUCT = Chem.MolFromSmarts('C=CC')
    MAP_NUMS = (InterMolecularReaction.NUCLEOPHILE_EAS_MAP_NUM, models.Template.EAS_MAP_NUM)

    def _get_template_mol(self, template):
        self.template = template.friedel_crafts_mol

    def _validate_reacting_atom(self, atom):
        if atom.GetSymbol() != 'C' or atom.GetTotalNumHs() == 0 or not atom.GetIsAromatic():
            raise InvalidMolecule(
                'The nucleophilic reacting atom must be an aromatic carbon with at least one hydrogen for a friedel crafts reaction!')

    def _validate_template(self):
        if self.template is None or not self.template.GetSubstructMatch(self.REQUIRED_SUBSTRUCT):
            raise InvalidMolecule(
                'The template molecule must contain a "C=CC" substructure for a friedel crafts reaction!')

    def _create_reactants(self):
        self._process_nucleophile()
        self.reactants = [self.reacting_mol, self.template]

    def _process_nucleophile(self):

        if isinstance(self.model, models.Sidechain):
            self._tag_sidechain_connection_atom()
        elif isinstance(self.model, models.Monomer):
            self._tag_monomer_connection_atom()
        else:
            raise InvalidMolecule(
                f'The nucleophile needs to be an instance of a sidechain or a monomer!')

    def _create_product(self):
        self.products = [utils.connect_mols(*self.reactants, map_nums=self.MAP_NUMS, clear_map_nums=False)]


class TsujiTrost(InterMolecularReaction):
    REQUIRED_SUBSTRUCT = Chem.MolFromSmarts('C=CC')
    MAP_NUMS = (InterMolecularReaction.NUCLEOPHILE_EAS_MAP_NUM, models.Template.EAS_MAP_NUM)

    def _get_template_mol(self, template):
        self.template = template.tsuji_trost_mol

    def _validate_reacting_atom(self, atom):
        if atom.GetSymbol() not in ('N', 'O', 'S') or atom.GetTotalNumHs() == 0:
            raise InvalidMolecule(
                'The nucleophilic reacting atom must be a heteroatom with at least one hydrogen for a tsuji trost reaction!')

    def _validate_template(self):
        if self.template is None or not self.template.GetSubstructMatch(self.REQUIRED_SUBSTRUCT):
            raise InvalidMolecule(
                'The template molecule must contain a "C=CC" substructure for a tsuji trost reaction!')

    def _create_reactants(self):
        self._process_nucleophile()
        self.reactants = [self.reacting_mol, self.template]

    def _process_nucleophile(self):

        if isinstance(self.model, models.Sidechain):
            self._tag_sidechain_connection_atom()
        elif isinstance(self.model, models.Monomer):
            self._tag_monomer_connection_atom()
        else:
            raise InvalidMolecule(
                f'The nucleophile needs to be an instance of a sidechain or a monomer!')

    def _create_product(self):
        self.products = [utils.connect_mols(*self.reactants, map_nums=self.MAP_NUMS, clear_map_nums=False)]


class PictetSpangler(InterMolecularReaction):
    REQUIRED_SUBSTRUCT = Chem.MolFromSmarts('C(=O)')
    BACKBONE_OLIGOMERIZATION_MAP_NUM = models.Backbone.MAP_NUM
    SIDECHAIN_OLIGOMERIZATION_MAP_NUM = InterMolecularReaction.NUCLEOPHILE_WC_MAP_NUM
    TEMPLATE_OLIGOMERIZATION_MAP_NUM = models.Template.OLIGOMERIZATION_MAP_NUM
    TEMPLATE_ALDEHYDE_O_MAP_NUM = models.Template.PS_OXYGEN_MAP_NUM
    TEMPLATE_ALDEHYDE_C_MAP_NUM = models.Template.PS_CARBON_MAP_NUM

    def _get_template_mol(self, template):
        self.template = template.pictet_spangler_mol

    def _validate_reacting_atom(self, atom):
        if atom.GetSymbol() != 'C' or atom.GetTotalNumHs() == 0 or not atom.GetIsAromatic():
            raise InvalidMolecule(
                'The nucleophilic reacting atom must be an aromatic carbon with at least one hydrogen for a pictet spangler reaction!')

        self.reacting_atom = atom

    def _validate_template(self):
        if self.template is None or not self.template.GetSubstructMatch(self.REQUIRED_SUBSTRUCT):
            raise InvalidMolecule(
                'The template molecule must contain a "C=CC" substructure for a tsuji trost reaction!')

    def _validate_sidechain(self):
        connection_atom = utils.find_atom(self.reacting_mol, self.SIDECHAIN_OLIGOMERIZATION_MAP_NUM)
        paths = Chem.FindAllPathsOfLengthN(self.reacting_mol, 3, useBonds=False,
                                           rootedAtAtom=self.reacting_atom.GetIdx())
        atoms = set().union([atom for path in paths for atom in path])
        if connection_atom.GetIdx() not in atoms - set(atom.GetIdx() for atom in connection_atom.GetNeighbors()):
            raise InvalidMolecule(
                'The attachment point of the reacting sidechain must be two atoms away from the reacting atom in a pictet spangler reaction!')

    def _validate_monomer(self):
        n_term_atom = utils.find_atom(self.reacting_mol, self.BACKBONE_NITROGEN_MAP_NUM)
        paths = Chem.FindAllPathsOfLengthN(self.reacting_mol, 5, useBonds=False,
                                           rootedAtAtom=self.reacting_atom.GetIdx())
        atoms = set().union([atom for path in paths for atom in path])
        if n_term_atom.GetIdx() not in atoms - set(atom.GetIdx() for atom in n_term_atom.GetNeighbors()):
            raise InvalidMolecule('The reacting atom in the monomer must be 4 atoms away from the N-terminus!')

    def _create_reactants(self):
        self._process_nucleophile()
        map_nums = (self.TEMPLATE_OLIGOMERIZATION_MAP_NUM, self.BACKBONE_NITROGEN_MAP_NUM)

        if isinstance(self.model, models.Sidechain):
            monomer = self._create_monomer()
            self.reactants = [utils.connect_mols(monomer, self.template, map_nums=map_nums, clear_map_nums=False)]
        else:
            # remove wildcard atom attached to n-terminus added in call to tag_monomer_connection_atom()
            wildcard_atom = utils.find_atom(self.reacting_atom, self.N_TERM_WILDCARD_MAP_NUM)
            self.reacting_mol = Chem.RWMol(self.reacting_mol)
            self.reacting_mol.RemoveAtom(wildcard_atom.GetIdx())

            # connect template to monomer
            self.reactants = [utils.connect_mols(self.reacting_mol, self.template,
                                                 map_nums=map_nums, clear_map_nums=False)]

    def _process_nucleophile(self):

        if isinstance(self.model, models.Sidechain):
            self._tag_sidechain_connection_atom(change_to_wildcard=False)
            self._validate_sidechain()
        elif isinstance(self.model, models.Monomer):
            self._tag_monomer_connection_atom()
            self._validate_monomer()
        else:
            raise InvalidMolecule(
                f'The nucleophile needs to be an instance of a sidechain or a monomer!')

    def _create_product(self):
        reactant = Chem.RWMol(self.reactants[0])  # copy of reactant molecule that will be transformed into product

        # reset atom map number on the INSTANCE VARIABLE reactant
        for atom in self.reactants[0].GetAtoms():
            if atom.GetAtomMapNum() == self.TEMPLATE_ALDEHYDE_O_MAP_NUM:
                atom.SetAtomMapNum(0)
                break

        # remove oxygen from unmasked aldehyde on LOCAL VARIABLE reactant
        for atom in list(reactant.GetAtoms()):
            if atom.GetAtomMapNum() == self.TEMPLATE_ALDEHYDE_O_MAP_NUM:
                reactant.RemoveAtom(atom.GetIdx())
                break

        # create bond between sidechain EAS carbon and unmasked aldehyde carbon
        map_nums = (self.NUCLEOPHILE_EAS_MAP_NUM, self.TEMPLATE_ALDEHYDE_C_MAP_NUM)
        reactant = utils.connect_mols(reactant, map_nums=map_nums, clear_map_nums=False)

        # create bond between peptide nitrogen and unmasked aldehyde carbon
        map_nums = (self.BACKBONE_NITROGEN_MAP_NUM, self.TEMPLATE_ALDEHYDE_C_MAP_NUM)
        self.products = [utils.connect_mols(reactant, map_nums=map_nums, clear_map_nums=False)]


class Pyrroloindoline(InterMolecularReaction):
    REQUIRED_SUBSTRUCT = Chem.MolFromSmarts('C=CC')
    INDOLE = Chem.MolFromSmarts('*c1c[nH]c2ccccc12')
    STEREOCHEMISTRIES = [('CCW', 'CCW'), ('CW', 'CW')]
    ADJ_CARBON_MAP_NUM = 7
    SIDECHAIN_OLIGOMERIZATION_MAP_NUM = InterMolecularReaction.NUCLEOPHILE_WC_MAP_NUM
    TEMPLATE_EAS_MAP_NUM = models.Template.EAS_MAP_NUM
    BACKBONE_OLIGOMERIZATION_MAP_NUM = models.Backbone.MAP_NUM

    def _get_template_mol(self, template):
        self.template = template.pyrroloindoline_mol

    def _validate_reacting_atom(self, atom):
        self._tag_sidechain_connection_atom(change_to_wildcard=False)

        num_neighbors_hydrogens = [neighbor.GetTotalNumHs() for neighbor in atom.GetNeighbors()]
        connection_atom = self._get_connection_atom()
        connection_atom_neighbors = [neighbor.GetIdx() for neighbor in connection_atom.GetNeighbors()]

        # reaction initiated at carbon that is the attachment point to amino acid backbone, which contains no hydrogens
        if atom.GetTotalNumHs() != 0 \
                or atom.GetSymbol() != 'C' \
                or 1 not in num_neighbors_hydrogens \
                or atom.GetIdx() not in connection_atom_neighbors:
            raise InvalidMolecule('The reacting atom of a pyrroloindoline reaction must be an aromatic carbon with no'
                                  ' hydrogens, adjacent to the connecting atom between the sidechain and backbone, and'
                                  ' with a neighbor that has exactly one hydrogen.')

    def _get_connection_atom(self):
        self._validate_reacting_mol()
        for atom in self.reacting_mol.GetAtoms():
            if atom.GetAtomMapNum() == self.SIDECHAIN_OLIGOMERIZATION_MAP_NUM:
                return atom

    def _validate_template(self):
        if self.template is None or not self.template.GetSubstructMatch(self.REQUIRED_SUBSTRUCT):
            raise InvalidMolecule(
                'The template molecule must contain a "C=CC" substructure for a pyrroloindoline reaction!')

    def _validate_reacting_mol(self):
        reacting_mol = Chem.MolFromSmiles(Chem.MolToSmiles(self.reacting_mol))
        if not reacting_mol.GetSubstructMatch(self.INDOLE) or not isinstance(self.model, models.Sidechain):
            raise InvalidMolecule(
                'The reacting molecule in a pyrroloindoline reaction must be a sidechain containing an indole!')

    def _create_reactants(self):
        monomer = self._create_monomer()

        # attach a methyl to n-terminus so reaction matches indole at any position in peptide chain
        extra_atom = Chem.MolFromSmarts(f'[CH4:{self.N_TERM_WILDCARD_MAP_NUM}]')
        map_nums = (self.BACKBONE_NITROGEN_MAP_NUM, self.N_TERM_WILDCARD_MAP_NUM)
        self.monomer = utils.connect_mols(monomer, extra_atom, map_nums=map_nums, clear_map_nums=False)

        for atom in self.monomer.GetAtoms():
            if atom.GetAtomMapNum() == self.N_TERM_WILDCARD_MAP_NUM:  # change methyl to wildcard
                utils.atom_to_wildcard(atom)
            elif atom.GetAtomMapNum() == self.NUCLEOPHILE_EAS_MAP_NUM:  # tag carbon adjacent to reacting atom
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetTotalNumHs() == 1 and neighbor.GetIsAromatic():
                        neighbor.SetAtomMapNum(self.ADJ_CARBON_MAP_NUM)

        self.reactants = [self.monomer, self.template]

    def _create_product(self):

        # copy of monomer that will be changed into product
        monomer = deepcopy(self.monomer)

        for atom in monomer.GetAtoms():
            if atom.GetAtomMapNum() == self.NUCLEOPHILE_EAS_MAP_NUM:
                eas_atom = atom
            elif atom.GetAtomMapNum() == self.ADJ_CARBON_MAP_NUM:
                adj_atom = atom

        # bond between eas_atom (i.e. the reacting atom) and adj_atom becomes single bond
        bond = monomer.GetBondBetweenAtoms(eas_atom.GetIdx(), adj_atom.GetIdx())
        bond.SetBondType(Chem.BondType.SINGLE)
        eas_atom.SetNumExplicitHs(eas_atom.GetTotalNumHs() + 1)
        adj_atom.SetNumExplicitHs(adj_atom.GetTotalNumHs() + 1)

        for stereo1, stereo2 in self.STEREOCHEMISTRIES:
            # merge backbone nitrogen to adjacent carbon
            map_nums = (self.BACKBONE_NITROGEN_MAP_NUM, self.ADJ_CARBON_MAP_NUM)
            reactant = utils.connect_mols(monomer, map_nums=map_nums, clear_map_nums=False, stereo=stereo1)

            # merge template with monomer
            map_nums = (self.TEMPLATE_EAS_MAP_NUM, self.NUCLEOPHILE_EAS_MAP_NUM)
            self.products.append(utils.connect_mols(reactant, self.template,
                                                    map_nums=map_nums, clear_map_nums=False, stereo=stereo2))


class IntraMolecularReaction:
    def generate(self, reacting_mol):
        self._initialize(reacting_mol)
        self._create_product()
        return create_reaction_smarts([self.reacting_mol], self.products)

    def _initialize(self, reacting_mol):
        self.reacting_mol = reacting_mol
        self.products = []

    def _create_product(self):
        pass


class TemplatePictetSpangler(IntraMolecularReaction):

    STEREOCHEMISTRY = 'CCW'
    ALDEHYDE_O_MAP_NUM = models.Template.PS_OXYGEN_MAP_NUM
    ALDEHYDE_C_MAP_NUM = models.Template.PS_CARBON_MAP_NUM
    EAS_MAP_NUM = models.Template.EAS_MAP_NUM
    NITROGEN_MAP_NUM = models.Template.TEMPLATE_PS_NITROGEN_MAP_NUM

    def _initialize(self, reacting_mol):
        super()._initialize(reacting_mol)

        try:
            self.reacting_mol = reacting_mol.template_pictet_spangler_mol
        except AttributeError:
            raise InvalidMolecule(
                'The reacting mol must be a template with a valid template pictet spangler mol property!')

    def _create_product(self):
        # copy of template that turns into product
        template = Chem.RWMol(self.reacting_mol)

        # reset atom map number on the INSTANCE VARIABLE template
        for atom in self.reacting_mol.GetAtoms():
            if atom.GetAtomMapNum() == self.ALDEHYDE_O_MAP_NUM:
                atom.SetAtomMapNum(0)
                break

        # remove oxygen from unmasked aldehyde on LOCAL VARIABLE template
        for atom in list(template.GetAtoms()):
            if atom.GetAtomMapNum() == self.ALDEHYDE_O_MAP_NUM:
                template.RemoveAtom(atom.GetIdx())
            elif atom.GetAtomMapNum() == self.ALDEHYDE_C_MAP_NUM:
                atom.SetNumExplicitHs(atom.GetTotalNumHs() + 2)

        # create bond between the peptide nitrogen atom and the aldehyde carbon
        map_nums = (self.NITROGEN_MAP_NUM, self.ALDEHYDE_C_MAP_NUM)
        template = utils.connect_mols(template, map_nums=map_nums, clear_map_nums=False)

        # create bond between reacting atom and the aldehyde carbon
        map_nums = (self.ALDEHYDE_C_MAP_NUM, self.EAS_MAP_NUM)
        self.products.append(utils.connect_mols(template, map_nums=map_nums,
                                                stereo=self.STEREOCHEMISTRY, clear_map_nums=False))
