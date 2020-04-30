from itertools import chain

from rdkit import Chem
from rdkit.Chem import AllChem

import cpmg.models as models
import macrocycles.utils as utils
from cpmg.exceptions import InvalidMolecule


class InterMolecularReaction:
    BACKBONE_NITROGEN_MAP_NUM = 100
    BACKBONE_CARBOXYL_MAP_NUM = 101
    N_TERM_WILDCARD_MAP_NUM = 102
    C_TERM_WILDCARD_MAP_NUM = 103

    def generate(self, reacting_mol, template, reacting_atom, model):
        self.initialize(reacting_mol, template, model)
        self._validate_reacting_atom(reacting_atom)
        self._validate_template()
        self._create_reactants()
        self._create_product()
        return self._create_reaction_smarts()

    def initialize(self, reacting_mol, template, model):
        self.model = model
        self.reacting_mol = reacting_mol
        self._get_template_mol(template)
        self.reactants = []
        self.product = None

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

    def _create_reaction_smarts(self):
        reactants = [Chem.MolToSmiles(reactant) for reactant in self.reactants]
        return '(' + '.'.join(reactants) + ')>>' + Chem.MolToSmiles(self.product)

    def _tag_sidechain_connection_atom(self, change_to_wildcard=True):

        matches = self.reacting_mol.GetSubstructMatches(Chem.MolFromSmarts('[CH3]'))
        for atom_idx in chain.from_iterable(matches):
            atom = self.reacting_mol.GetAtomWithIdx(atom_idx)
            if atom.GetIsotope() == 13:  # isotope labeling for methyl carbons that arent candidate attachments
                atom.SetIsotope(0)
            else:
                atom.SetAtomMapNum(self.NUCLEOPHILE_WC_MAP_NUM)  # pylint: disable=maybe-no-member
                if change_to_wildcard:
                    utils.atom_to_wildcard(atom)

    def _tag_monomer_connection_atom(self):

        n_term_patt1 = Chem.MolFromSmarts('[NH2]')  # regular amino acid backbone type
        n_term_patt2 = Chem.MolFromSmarts('[NH;R]')  # proline backbone type
        c_term_patt = Chem.MolFromSmarts('C(=O)[OH]')
        n_term_replacement1 = Chem.MolFromSmarts(
            f'[NH1:{self.BACKBONE_NITROGEN_MAP_NUM}][*:{self.N_TERM_WILDCARD_MAP_NUM}]')
        n_term_replacement2 = Chem.MolFromSmarts(
            f'[N;R{self.BACKBONE_NITROGEN_MAP_NUM}][*:{self.N_TERM_WILDCARD_MAP_NUM}]')
        c_term_replacement = Chem.MolFromSmarts(
            f'[C:{self.BACKBONE_CARBOXYL_MAP_NUM}](=O)[*:{self.C_TERM_WILDCARD_MAP_NUM}]')

        # if nucleophile doesnt have the specified pattern, ReplaceSubstracts() returns the mol unmodified
        self.reacting_mol = AllChem.ReplaceSubstructs(self.reacting_mol, n_term_patt1, n_term_replacement1)[0]
        self.reacting_mol = AllChem.ReplaceSubstructs(self.reacting_mol, n_term_patt2, n_term_replacement2)[0]
        self.reacting_mol = AllChem.ReplaceSubstructs(self.reacting_mol, c_term_patt, c_term_replacement)[0]


class FriedelCrafts(InterMolecularReaction):
    REQUIRED_SUBSTRUCT = Chem.MolFromSmarts('C=CC')
    NUCLEOPHILE_EAS_MAP_NUM = 1
    NUCLEOPHILE_WC_MAP_NUM = 2
    TEMPLATE_EAS_MAP_NUM = models.Template.FRIEDEL_CRAFTS_EAS_MAP_NUM
    MAP_NUMS = (NUCLEOPHILE_EAS_MAP_NUM, TEMPLATE_EAS_MAP_NUM)

    def _get_template_mol(self, template):
        self.template = template.friedel_crafts_mol

    def _validate_reacting_atom(self, atom):
        if atom.GetSymbol() != 'C' \
                or atom.GetTotalNumHs() == 0 \
                or not atom.GetIsAromatic():
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
        self.product = utils.connect_mols(*self.reactants, map_nums=self.MAP_NUMS, clear_map_nums=False)
