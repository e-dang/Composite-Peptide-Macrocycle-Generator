
import multiprocessing
from collections import namedtuple
from copy import deepcopy
from itertools import chain, product
from logging import INFO

from rdkit import Chem
from rdkit.Chem import AllChem

import macrocycles.config as config
import macrocycles.utils as utils

LOGGER = utils.create_logger(__name__, INFO)


########################################################################################################################
########################################################################################################################
########################################################################################################################

ReactionInfo = namedtuple('ReactionInfo', 'binary rxn_atom_idx type')


class Reaction():

    def __init__(self, validation, product_generator, name, reacting_atom, *reactants):

        self.validation = validation
        self.product_generator = product_generator
        self.name = name
        self.reacting_atom = reacting_atom
        self.reactants = reactants
        self.valid = False

    def __repr__(self):
        name = f'name: {self.name}'
        reactants = f'reactants: {[Chem.MolToSmiles(reactant) for reactant in self.reactants]}'
        prod = f'product: {Chem.MolToSmiles(self.product)}'
        return '\n'.join([name, reactants, prod])

    def __str__(self):
        if self.product is None:
            return None

        reactants = [Chem.MolToSmiles(reactant) for reactant in self.reactants]
        return '(' + '.'.join(reactants) + ')>>' + Chem.MolToSmiles(self.product)

    def __bool__(self):
        return self.validation()

    @property
    def product(self):
        try:
            return self.cached_product  # pylint: disable=access-member-before-definition
        except AttributeError:
            pass

        if self:
            self.cached_product = self.product_generator()
            return self.cached_product

        self.cached_product = None
        return self.cached_product

    @property
    def reaction(self):
        if self.product is None:
            return None

        return AllChem.ReactionFromSmarts(str(self))

    @property
    def binary(self):
        return self.reaction.ToBinary()

    @property
    def data(self):
        if self.product is None:
            return None

        return ReactionInfo(self.binary, self.reacting_atom.GetIdx(), self.name)

########################################################################################################################
########################################################################################################################
########################################################################################################################


class FriedelCafts(Reaction):

    def __init__(self, side_chain, template, reacting_atom):

        super().__init__(self.is_valid, self.generate_product, 'friedel_crafts', reacting_atom, side_chain, template)

    def is_valid(self):

        if self.valid:
            return True

        if self.reacting_atom.GetSymbol() == 'C' \
                and self.reacting_atom.GetTotalNumHs() != 0 \
                and self.reacting_atom.GetIsAromatic():
            self.preprocess()
            return True

        return False

    def preprocess(self):
        self.reactants = [self.reactants[0], Chem.MolFromSmarts(
            f'[*:{config.TEMP_WILDCARD_MAP_NUM}]/C=C/[CH3:{config.TEMP_EAS_MAP_NUM}]')]

    def generate_product(self):
        ignored_map_nums = [config.TEMP_WILDCARD_MAP_NUM, config.SC_WILDCARD_MAP_NUM]
        return utils.Base.merge(*self.reactants, ignored_map_nums=ignored_map_nums, clear_map_nums=False)

########################################################################################################################
########################################################################################################################
########################################################################################################################


class TsujiTrost(Reaction):

    def __init__(self, side_chain, template, reacting_atom):

        super().__init__(self.is_valid, self.generate_product, 'tsuji_trost', reacting_atom, side_chain, template)

    def is_valid(self):
        if self.valid:
            return True

        if self.reacting_atom.GetSymbol() in ['N', 'O', 'S'] and self.reacting_atom.GetTotalNumHs() > 0:
            self.preprocess()
            return True

        return False

    def preprocess(self):
        self.reactants = [self.reactants[0], Chem.MolFromSmarts(
            f'[*:{config.TEMP_WILDCARD_MAP_NUM}]/C=C/[CH3:{config.TEMP_EAS_MAP_NUM}]')]

    def generate_product(self):
        ignored_map_nums = [config.TEMP_WILDCARD_MAP_NUM, config.SC_WILDCARD_MAP_NUM]
        return utils.Base.merge(*self.reactants, ignored_map_nums=ignored_map_nums, clear_map_nums=False)

########################################################################################################################
########################################################################################################################
########################################################################################################################


class PictetSpangler(Reaction):
    CARBON_PS_MAP_NUM = 5
    CARBON_PEP_MAP_NUM = 6
    NITROGEN_PS_MAP_NUM = 7
    OXYGEN_PS_MAP_NUM = 8
    CARBOX_WILDCARD_MAP_NUM = 9
    ALKYN_WILDCARD_MAP_NUM = 10

    def __init__(self, side_chain, template, reacting_atom):

        super().__init__(self.is_valid, self.generate_product, 'pictet_spangler', reacting_atom, side_chain, template)

    def is_valid(self):
        if self.valid:
            return True

        side_chain, template = self.reactants
        wild_card = [atom for atom in side_chain.GetAtoms() if atom.GetAtomicNum() == 0][0]
        paths = Chem.FindAllPathsOfLengthN(side_chain, 3, useBonds=False, rootedAtAtom=self.reacting_atom.GetIdx())
        atoms = set().union([atom for path in paths for atom in path])

        # check reacting atom is aromatic carbon
        if not self.reacting_atom.GetIsAromatic() \
                or not self.reacting_atom.GetSymbol() == 'C' \
                or self.reacting_atom.GetTotalNumHs() == 0:
            return False

        # check that attachment point of side_chain is two atoms away from reacting atom
        if wild_card.GetIdx() not in atoms - set(atom.GetIdx() for atom in wild_card.GetNeighbors()):
            return False

        # check that template has unmasked aldehyde (this aldehyde will participate in pictet spangler)
        match = template.GetSubstructMatch(Chem.MolFromSmarts('C(=O)'))
        if not match:
            return False

        self.preprocess(wild_card, match)
        self.valid = True
        return True

    def preprocess(self, wild_card, substruct_match):
        monomers = self.create_monomers(wild_card)
        template = self.modify_template(substruct_match)
        self.create_reactant(monomers, template)

    def create_monomers(self, wild_card):

        # change side_chain wild_card back to carbon
        wild_card.SetAtomicNum(6)

        # get backbone and side_chain
        db = utils.MongoDataBase(logger=None)
        side_chain = self.reactants[0]
        backbone = Chem.Mol(db['molecules'].find_one({'_id': 'alpha'})['binary'])

        # tag backbone nitrogen and remove carboxyl group
        backbone = AllChem.ReplaceSubstructs(backbone, Chem.MolFromSmarts(
            'C(=O)O'), Chem.MolFromSmarts(f'[*:{self.CARBOX_WILDCARD_MAP_NUM}]'))[0]
        for atom in chain.from_iterable(backbone.GetSubstructMatches(Chem.MolFromSmarts('[NH2]'))):
            atom = backbone.GetAtomWithIdx(atom)
            if atom.GetSymbol() == 'N':
                atom.SetAtomMapNum(self.NITROGEN_PS_MAP_NUM)
                break

        # create monomer
        ignored_map_nums = [config.SC_EAS_MAP_NUM, self.NITROGEN_PS_MAP_NUM, self.CARBOX_WILDCARD_MAP_NUM]
        return utils.Base.merge(side_chain, backbone, ignored_map_nums=ignored_map_nums)

    def modify_template(self, substruct_match):
        template = self.reactants[1]

        # tag substruct match in template
        for atom in substruct_match:
            atom = template.GetAtomWithIdx(atom)
            if atom.GetSymbol() == 'C':
                atom.SetAtomMapNum(self.CARBON_PS_MAP_NUM)
            elif atom.GetSymbol() == 'O':
                atom.SetAtomMapNum(self.OXYGEN_PS_MAP_NUM)

        # change wild card atom to carbonyl
        unique = set()
        rxn = AllChem.ReactionFromSmarts('[#0:1][*:2]>>[C:1](=O)[*:2]')
        for prod in chain.from_iterable(rxn.RunReactants((template,))):
            Chem.SanitizeMol(prod)
            unique.add(prod.ToBinary())

        # set atom map number of peptide carbonyl
        template = Chem.Mol(unique.pop())
        for atom in chain.from_iterable(template.GetSubstructMatches(Chem.MolFromSmarts('C(=O)'))):
            atom = template.GetAtomWithIdx(atom)
            if atom.GetSymbol() == 'C' and atom.GetAtomMapNum() != self.CARBON_PS_MAP_NUM:
                atom.SetAtomMapNum(self.CARBON_PEP_MAP_NUM)

        # change cinnamoyl unit to wildcard
        unique = set()
        rxn = AllChem.ReactionFromSmarts(
            f'C/C=C/c1cc[$(c(F)),$(c)]c([*:{config.TEMP_EAS_MAP_NUM}])c1>>*[*:{config.TEMP_EAS_MAP_NUM}]')
        for prod in chain.from_iterable(rxn.RunReactants((template,))):
            Chem.SanitizeMol(prod)
            unique.add(prod.ToBinary())

        # tag new wildcard atom
        template = Chem.Mol(unique.pop())
        for atom in template.GetAtoms():
            if atom.GetAtomicNum() == 0:
                atom.SetAtomMapNum(config.TEMP_EAS_MAP_NUM)
                break

        # change alkyne to anther wildcard
        wild_card2 = Chem.MolFromSmarts(f'[*:{self.ALKYN_WILDCARD_MAP_NUM}]')
        alkyne = Chem.MolFromSmarts('C#C')
        template = AllChem.ReplaceSubstructs(template, alkyne, wild_card2)[0]

        return template

    def create_reactant(self, monomer, template):

        ignored_map_nums = [config.TEMP_EAS_MAP_NUM, config.SC_EAS_MAP_NUM, self.CARBOX_WILDCARD_MAP_NUM,
                            self.CARBON_PS_MAP_NUM, self.OXYGEN_PS_MAP_NUM, self.ALKYN_WILDCARD_MAP_NUM]
        self.reactants = [utils.Base.merge(monomer, template, ignored_map_nums=ignored_map_nums, clear_map_nums=False)]

    def generate_product(self):
        reactant = Chem.RWMol(self.reactants[0])

        # reset atom map number on reactant
        for atom in self.reactants[0].GetAtoms():
            if atom.GetAtomMapNum() == self.OXYGEN_PS_MAP_NUM:
                atom.SetAtomMapNum(0)

        # remove oxygen from unmasked aldehyde
        for atom in list(reactant.GetAtoms()):
            if atom.GetAtomMapNum() == self.OXYGEN_PS_MAP_NUM:
                reactant.RemoveAtom(atom.GetIdx())
                break

        # create bond between side_chain EAS carbon and unmasked aldehyde carbon
        ignored_map_nums = [config.TEMP_EAS_MAP_NUM, self.CARBON_PEP_MAP_NUM, self.ALKYN_WILDCARD_MAP_NUM,
                            self.NITROGEN_PS_MAP_NUM, self.CARBOX_WILDCARD_MAP_NUM]
        reactant = utils.Base.merge(reactant, ignored_map_nums=ignored_map_nums, clear_map_nums=False)

        # create bond between peptide nitrogen and unmasked aldehyde carbon
        ignored_map_nums = [config.TEMP_EAS_MAP_NUM, config.SC_EAS_MAP_NUM, self.ALKYN_WILDCARD_MAP_NUM,
                            self.CARBON_PEP_MAP_NUM, self.CARBOX_WILDCARD_MAP_NUM]
        return utils.Base.merge(reactant, ignored_map_nums=ignored_map_nums, clear_map_nums=False)


########################################################################################################################
########################################################################################################################
########################################################################################################################


NewReactions = namedtuple('NewReactions', 'reactions side_chain template')


class ReactionGenerator(utils.Base):
    """
    Class for generating atom mapped reaction SMARTS strings from atom mapped template and side chain SMILES strings.
    Inherits from Base.

    Attributes:
        side_chains (list): Contains the different atom mapped side chains.
        templates (list): Contains the different atom mapped templates.
        reaction (str): The type of reaction that is being performed (i.e. Friedel-Crafts, Pictet-Spengler,..).
    """

    _defaults = config.DEFAULTS['ReactionGenerator']

    def __init__(self, logger=LOGGER, make_db_connection=True):
        """
        Initializer.

        Args:
        """

        # I/O
        super().__init__(LOGGER, make_db_connection)

        # data
        self.side_chains = []
        self.templates = []

    def save_data(self):

        params = self._defaults['outputs']

        return self.to_mongo(params['col_reactions'])

    def load_data(self):
        """
        Overloaded method for loading input data.

        Returns:
            bool: True if successful.
        """

        try:
            params = self._defaults['inputs']
            self.side_chains = list(self.from_mongo(params['col_side_chains'], {
                'type': 'side_chain', 'connection': 'methyl'}))
            self.templates = list(self.from_mongo(params['col_templates'], {'type': 'template'}))
        except Exception:
            self.logger.exception('Unexpected exception occured.')
        else:
            return True

        return False

    def generate(self, reactions=[FriedelCafts, TsujiTrost, PictetSpangler]):

        try:
            dependent_rxns, independent_rxns = self.classify_reactions(reactions)
            dependent_args = product(self.side_chains, dependent_rxns, self.templates)
            independent_args = product(self.side_chains, independent_rxns)

            with multiprocessing.Pool() as pool:
                results = pool.starmap_async(ReactionGenerator.create_reaction, dependent_args)
                for rxns, side_chain, template in results.get():
                    self.accumulate_data(rxns, side_chain, template)

                results = pool.starmap_async(ReactionGenerator.create_reaction, independent_args)
                for rxns, side_chain, template in results.get():
                    self.accumulate_data(rxns, side_chain, template)
        except Exception:
            raise
        else:
            return True

        return False

    def generate_serial(self, reactions=[FriedelCafts, TsujiTrost, PictetSpangler]):

        try:
            dependent_rxns, independent_rxns = self.classify_reactions(reactions)

            for side_chain in self.side_chains:
                for template in self.templates:
                    for reaction in dependent_rxns:
                        rxn, _, _ = ReactionGenerator.create_reaction(side_chain, reaction, template)
                        self.accumulate_data(rxn, side_chain, template)
                for reaction in independent_rxns:
                    rxn, _, template = ReactionGenerator.create_reaction(side_chain, reaction)
                    self.accumulate_data(rxn, side_chain, template)
        except Exception:
            raise
        else:
            return True

        return False

    def generate_from_ids(self, side_chain_ids, template_ids, reactions=[FriedelCafts, TsujiTrost, PictetSpangler]):

        try:
            params = self._defaults['inputs']
            self.side_chains = list(self.from_mongo(params['col_side_chains'], {
                'type': 'side_chain', '_id': {'$in': side_chain_ids}}))
            self.templates = list(self.from_mongo(params['col_templates'], {
                'type': 'template', '_id': {'$in': template_ids}}))
            dependent_rxns, independent_rxns = self.classify_reactions(reactions)

            for side_chain in self.side_chains:
                for template in self.templates:
                    for reaction in dependent_rxns:
                        rxn, _, _ = ReactionGenerator.create_reaction(side_chain, reaction, template)
                        self.accumulate_data(rxn, side_chain, template)

                for reaction in independent_rxns:
                    rxn, _, template = ReactionGenerator.create_reaction(side_chain, reaction)
                    self.accumulate_data(rxn, side_chain, template)
        except Exception:
            raise
        else:
            return True

        return False

    def accumulate_data(self, reactions, side_chain, template):
        """
        Stores all data associated with the different reactions into a dictionary and appends it so self.result_data.

        Args:
            template (dict): The associated data of the templates.
            side_chain (dict): The associated data of the side chain that the regioisomers are derived from.
            reaction (str): The atom mapped reaction SMARTS string.
        """

        if template is not None:
            chunk = len(reactions) * self.templates.index(template)
            applicable_template = template['_id']
        else:
            chunk = len(reactions)
            applicable_template = 'all'

        for i, (smarts, (binary, rxn_atom_idx, rxn_type)) in enumerate(reactions.items()):
            doc = {'_id': side_chain['_id'] + str(chunk + i) + rxn_type[0],
                   'type': rxn_type,
                   'binary': binary,
                   'smarts': smarts,
                   'rxn_atom_idx': rxn_atom_idx,
                   'template': applicable_template,
                   'side_chain': {'_id': side_chain['_id'],
                                  'parent_side_chain': side_chain['parent_side_chain']['_id'],
                                  'conn_atom_idx': side_chain['conn_atom_idx']}}
            self.result_data.append(doc)

    @staticmethod
    def create_reaction(side_chain, reaction, template=None):

        try:
            sc_mol = Chem.Mol(side_chain['binary'])
            temp_mol = Chem.Mol(template['binary'])
        except TypeError:
            temp_mol = None

        reactions = {}
        ReactionGenerator.tag_wildcard_atom(sc_mol)

        for atom in sc_mol.GetAtoms():

            if atom.GetAtomMapNum() == config.SC_WILDCARD_MAP_NUM:
                continue

            atom.SetAtomMapNum(config.SC_EAS_MAP_NUM)
            rxn = reaction(deepcopy(sc_mol), deepcopy(temp_mol), atom)
            if rxn:
                reactions[str(rxn)] = rxn.data

            atom.SetAtomMapNum(0)

        return NewReactions(reactions, side_chain, template)

    @staticmethod
    def tag_wildcard_atom(side_chain):
        matches = side_chain.GetSubstructMatches(Chem.MolFromSmarts('[CH3]'))
        for pair in matches:
            for atom_idx in pair:
                atom = side_chain.GetAtomWithIdx(atom_idx)
                atom.SetAtomMapNum(config.SC_WILDCARD_MAP_NUM)
                utils.atom_to_wildcard(atom)
                return atom

    def classify_reactions(self, reactions):
        params = self._defaults['template_dependent_rxns']
        dependent_rxns, independent_rxns = [], []
        for func in reactions:
            dependent_rxns.append(func) if func.__name__ in params else independent_rxns.append(func)

        return dependent_rxns, independent_rxns

    # def reaction_from_vector(mol, vector=None):

    #     mol = mol.split('*')
    #     for i, substr in enumerate(mol):
    #         if i == 0:
    #             mol = substr
    #         else:
    #             mol = f'[*:{i}]'.join([mol, substr])

    #     print(mol)
    #     mol = Chem.MolFromSmiles(mol)

    #     idxs = [idx for idx_tup in vector if 0 not in idx_tup for idx in idx_tup]
    #     idxs = set(idxs)
    #     print(idxs)
    #     for i, idx in enumerate(idxs, start=i + 1):
    #         mol.GetAtomWithIdx(idx).SetAtomMapNum(i)

    #     old_mol = deepcopy(mol)
    #     mol = Chem.RWMol(mol)
    #     idxs = set()
    #     for idx_tup in vector:
    #         idx1, idx2 = idx_tup
    #         idxs.add(idx1)
    #         idxs.add(idx2)

    #         if 0 in idx_tup:
    #             idx = idx_tup[0] if idx_tup[0] != 0 else idx_tup[1]
    #             mol.RemoveAtom(idx)
    #         else:
    #             idx1, idx2 = idx_tup
    #             idxs.add(idx1)
    #             idxs.add(idx2)
    #             mol.AddBond(idx1, idx2, Chem.rdchem.BondType.SINGLE)

    #     reaction = '(' + Chem.MolToSmiles(old_mol) + ')>>' + Chem.MolToSmiles(mol)

    #     return reaction
