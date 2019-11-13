
import multiprocessing
from collections import namedtuple
from copy import deepcopy
from itertools import chain, product
from logging import INFO

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

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

    def sc_attach_adjacent_N(self, carboxyl_map_num, nitrogen_map_num, clear_map_nums=True):

        # get side_chain and backbone
        side_chain = self.reactants[0]
        db = utils.MongoDataBase(logger=None)
        backbone = Chem.Mol(db['molecules'].find_one({'_id': 'alpha'})['binary'])

        # tag backbone nitrogen and remove carboxyl group
        backbone = AllChem.ReplaceSubstructs(backbone, Chem.MolFromSmarts(
            'C(=O)O'), Chem.MolFromSmarts(f'[*:{carboxyl_map_num}]'))[0]
        for atom in chain.from_iterable(backbone.GetSubstructMatches(Chem.MolFromSmarts('[NH2]'))):
            atom = backbone.GetAtomWithIdx(atom)
            if atom.GetSymbol() == 'N':
                atom.SetAtomMapNum(nitrogen_map_num)
                break

        # create monomer with side_chain adjacent to nitrogen
        ignored_map_nums = [config.SC_EAS_MAP_NUM, nitrogen_map_num, carboxyl_map_num]
        return utils.Base.merge(side_chain, backbone, ignored_map_nums=ignored_map_nums, clear_map_nums=clear_map_nums)

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
    C_TERM_WILDCARD_MAP_NUM = 9
    ALKYN_WILDCARD_MAP_NUM = 10

    def __init__(self, side_chain, template, reacting_atom):

        super().__init__(self.is_valid, self.generate_product, 'pictet_spangler', reacting_atom, side_chain, template)

    def is_valid(self):
        if self.valid:
            return True

        side_chain, template = self.reactants

        # check reacting atom is aromatic carbon
        if not self.reacting_atom.GetIsAromatic() \
                or self.reacting_atom.GetSymbol() != 'C' \
                or self.reacting_atom.GetTotalNumHs() == 0:
            return False

        # check that attachment point of side_chain is two atoms away from reacting atom
        wild_card = [atom for atom in side_chain.GetAtoms() if atom.GetAtomicNum() == 0][0]
        paths = Chem.FindAllPathsOfLengthN(side_chain, 3, useBonds=False, rootedAtAtom=self.reacting_atom.GetIdx())
        atoms = set().union([atom for path in paths for atom in path])
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
        monomers = self.modify_side_chain(wild_card)
        template = self.modify_template(substruct_match)
        self.create_reactant(monomers, template)

    def modify_side_chain(self, wild_card):

        # change side_chain wild_card back to carbon
        wild_card.SetAtomicNum(6)

        # create monomer
        return super().sc_attach_adjacent_N(self.C_TERM_WILDCARD_MAP_NUM, self.NITROGEN_PS_MAP_NUM)

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

        ignored_map_nums = [config.TEMP_EAS_MAP_NUM, config.SC_EAS_MAP_NUM, self.C_TERM_WILDCARD_MAP_NUM,
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
                            self.NITROGEN_PS_MAP_NUM, self.C_TERM_WILDCARD_MAP_NUM]
        reactant = utils.Base.merge(reactant, ignored_map_nums=ignored_map_nums, clear_map_nums=False)

        # create bond between peptide nitrogen and unmasked aldehyde carbon
        ignored_map_nums = [config.TEMP_EAS_MAP_NUM, config.SC_EAS_MAP_NUM, self.ALKYN_WILDCARD_MAP_NUM,
                            self.CARBON_PEP_MAP_NUM, self.C_TERM_WILDCARD_MAP_NUM]
        return utils.Base.merge(reactant, ignored_map_nums=ignored_map_nums, clear_map_nums=False)

########################################################################################################################
########################################################################################################################
########################################################################################################################


class PyrroloIndolene(Reaction):

    ADJ_CARBON_MAP_NUM = 5
    NITROGEN_MAP_NUM = 6
    C_TERM_WILDCARD_MAP_NUM = 7
    N_TERM_WILDCARD_MAP_NUM = 8

    def __init__(self, side_chain, template, reacting_atom):

        super().__init__(self.is_valid, self.generate_product, 'pyrrolo_indolene', reacting_atom, side_chain, template)

    def is_valid(self):
        if self.valid:
            return True

        side_chain, template = self.reactants

        # side_chain must be an indole
        match = side_chain.GetSubstructMatch(Chem.MolFromSmarts('*c1c[nH]c2ccccc12'))
        if not match:
            return False

        # reaction initiated at carbon that is the attachment point to amino acid backbone, which contains no hydrogens
        if self.reacting_atom.GetTotalNumHs() != 0 or self.reacting_atom.GetSymbol() != 'C':
            return False

        # reaction also involves carbon adjacent to reacting_atom which must have a hydrogen
        if 1 not in [atom.GetTotalNumHs() for atom in self.reacting_atom.GetNeighbors()]:
            return False

        # check if carbon is the correct bridged carbon (adjacent to nitrogen)
        wild_card = [atom for atom in side_chain.GetAtoms() if atom.GetAtomicNum() == 0][0]
        if self.reacting_atom.GetIdx() not in [atom.GetIdx() for atom in wild_card.GetNeighbors()]:
            return False

        self.preprocess()
        self.valid = True
        return True

    def preprocess(self):
        side_chain = self.modify_side_chain()
        self.reactants = [side_chain, Chem.MolFromSmarts(
            f'[*:{config.TEMP_WILDCARD_MAP_NUM}]/C=C/[CH3:{config.TEMP_EAS_MAP_NUM}]')]

    def modify_side_chain(self):

        # change side_chain wild card back to carbon
        for atom in self.reactants[0].GetAtoms():
            if atom.GetAtomicNum() == 0:
                atom.SetAtomicNum(6)
                break

        # get side_chain and backbone
        side_chain = super().sc_attach_adjacent_N(self.C_TERM_WILDCARD_MAP_NUM, self.NITROGEN_MAP_NUM)

        # add wildcard atom to backbone nitrogen
        unique = set()
        rxn = AllChem.ReactionFromSmarts('[NH2:1][*:2]>>*[NH1:1][*:2]')
        for prod in chain.from_iterable(rxn.RunReactants((side_chain,))):
            Chem.SanitizeMol(prod)
            unique.add(prod.ToBinary())

        # tag new wildcard atom, peptide nitrogen, and carbon adjacent to reacting atom that has a hydrogen
        side_chain = Chem.Mol(unique.pop())
        for atom in side_chain.GetAtoms():
            if atom.GetAtomMapNum() == 0 and atom.GetAtomicNum() == 0:
                atom.SetAtomMapNum(self.N_TERM_WILDCARD_MAP_NUM)
                neighbor = atom.GetNeighbors()[0]
                neighbor.SetAtomMapNum(self.NITROGEN_MAP_NUM)
            elif atom.GetAtomMapNum() == config.SC_EAS_MAP_NUM:
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetTotalNumHs() == 1 and neighbor.GetIsAromatic():
                        neighbor.SetAtomMapNum(self.ADJ_CARBON_MAP_NUM)
                        break

        return side_chain

    def generate_product(self):

        side_chain, template = deepcopy(self.reactants)

        atoms = []
        for atom in side_chain.GetAtoms():
            if atom.GetAtomMapNum() in (config.SC_EAS_MAP_NUM, self.ADJ_CARBON_MAP_NUM):
                atoms.append(atom)
                if len(atoms) == 2:
                    break

        # bond between reacting_atom and adjacent atom becomes single bond
        bond = side_chain.GetBondBetweenAtoms(atoms[0].GetIdx(), atoms[1].GetIdx())
        bond.SetBondType(Chem.BondType.SINGLE)
        atoms[0].SetNumExplicitHs(atoms[0].GetTotalNumHs() + 1)
        atoms[1].SetNumExplicitHs(atoms[1].GetTotalNumHs() + 1)

        # merge peptide nitrogen to adj carbon
        ignored_map_nums = [self.C_TERM_WILDCARD_MAP_NUM, config.SC_EAS_MAP_NUM, self.N_TERM_WILDCARD_MAP_NUM]
        side_chain = utils.Base.merge(side_chain, ignored_map_nums=ignored_map_nums, clear_map_nums=False)

        # merge template with side_chain
        ignored_map_nums = [self.C_TERM_WILDCARD_MAP_NUM, self.N_TERM_WILDCARD_MAP_NUM,
                            config.TEMP_WILDCARD_MAP_NUM, self.ADJ_CARBON_MAP_NUM, self.NITROGEN_MAP_NUM]
        return utils.Base.merge(side_chain, template, ignored_map_nums=ignored_map_nums, clear_map_nums=False)

########################################################################################################################
########################################################################################################################
########################################################################################################################


class TemplateOnlyPictetSpangler(Reaction):

    NITROGEN_PS_MAP_NUM = 5
    CARBON_PS_MAP_NUM = 6
    ALDEHYDE_CARBON = 7
    ALDEHYDE_OXYGEN = 8

    def __init__(self, template):
        super().__init__(self.is_valid, self.generate_product, 'template_only_pictet_spangler', None, template)

    def is_valid(self):
        if self.valid:  # pylint: disable=access-member-before-definition
            return True

        self.find_reacting_atom()

        # a valid reacting atom could not be found
        if self.reacting_atom is None:
            return False

        self.preprocess()
        self.valid = True
        return True

    def find_reacting_atom(self):

        template = self.reactants[0]

        # template must contain unmasked aldehyde
        for atom in chain.from_iterable(template.GetSubstructMatches(Chem.MolFromSmarts('C(=O)'))):
            atom = template.GetAtomWithIdx(atom)
            if atom.GetSymbol() == 'C':
                atom.SetAtomMapNum(self.ALDEHYDE_CARBON)
                aldehyde_carbon = atom.GetIdx()
            else:
                atom.SetAtomMapNum(self.ALDEHYDE_OXYGEN)
        try:
            # find all atoms 5 atoms away from aldehyde_carbon
            paths = Chem.FindAllPathsOfLengthN(template, 5, useBonds=False, rootedAtAtom=aldehyde_carbon)
            candidate_atoms = set().union([template.GetAtomWithIdx(atom) for path in paths for atom in path])

            # reacting atom must be on phenyl ring ortho to the unmasked aldehyde and para to the alkene where EAS takes
            # place, as well as have one hydrogen
            candidate_atoms = filter(lambda x: x.GetIsAromatic() and x.GetTotalNumHs() == 1, candidate_atoms)
            for atom in candidate_atoms:
                neighbor_hydrogen_counts = 0

                for neighbor in atom.GetNeighbors():
                    neighbor_hydrogen_counts += neighbor.GetTotalNumHs()

                if neighbor_hydrogen_counts == 1:
                    atom.SetAtomMapNum(self.CARBON_PS_MAP_NUM)
                    self.reacting_atom = atom
                    break
        except UnboundLocalError:  # no aldehyde present in template
            pass

    def preprocess(self):

        template = self.reactants[0]  # pylint: disable=access-member-before-definition

        # change wild card atom to amide
        unique = set()
        rxn = AllChem.ReactionFromSmarts('[#0:1][*:2]>>[*:2][C:1](=O)[NH1]*')
        for prod in chain.from_iterable(rxn.RunReactants((template,))):
            Chem.SanitizeMol(prod)
            unique.add(prod.ToBinary())

        # tag amide nitrogen and wild card atoms
        template = Chem.Mol(unique.pop())
        for atom in template.GetAtoms():
            if atom.GetSymbol() == 'N':
                atom.SetAtomMapNum(self.NITROGEN_PS_MAP_NUM)
            elif atom.GetAtomicNum() == 0 and atom.GetAtomMapNum() == 0:
                atom.SetAtomMapNum(config.TEMP_WILDCARD_MAP_NUM)

        # change alkene to wild card
        unique = set()
        rxn = AllChem.ReactionFromSmarts('[*:2]C=CC>>[*:2]*')
        for prod in chain.from_iterable(rxn.RunReactants((template,))):
            Chem.SanitizeMol(prod)
            unique.add(prod.ToBinary())

        # tag new wild card atom
        template = Chem.Mol(unique.pop())
        for atom in template.GetAtoms():
            if atom.GetAtomMapNum() == 0 and atom.GetAtomicNum() == 0:
                atom.SetAtomMapNum(config.TEMP_EAS_MAP_NUM)

        self.reactants = [template]

    def generate_product(self):

        template = Chem.RWMol(self.reactants[0])

        # remove oxygen atom map number on reactant
        for atom in self.reactants[0].GetAtoms():
            if atom.GetAtomMapNum() == self.ALDEHYDE_OXYGEN:
                atom.SetAtomMapNum(0)
                break

        # remove aldehyde oxygen atom
        for atom in list(template.GetAtoms()):
            if atom.GetAtomMapNum() == self.ALDEHYDE_OXYGEN:
                template.RemoveAtom(atom.GetIdx())
                break

        # create bond between template pictet_spangler carbon and aldehyde carbon
        ignored_map_nums = [config.TEMP_WILDCARD_MAP_NUM, config.TEMP_EAS_MAP_NUM, self.NITROGEN_PS_MAP_NUM]
        template = utils.Base.merge(template, ignored_map_nums=ignored_map_nums, clear_map_nums=False)

        # create bond between peptide nitrogen and unmasked aldehyde carbon
        ignored_map_nums = [config.TEMP_WILDCARD_MAP_NUM, config.TEMP_EAS_MAP_NUM, self.CARBON_PS_MAP_NUM]
        return utils.Base.merge(template, ignored_map_nums=ignored_map_nums, clear_map_nums=False)


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

    def save_data(self, to_mongo=True, to_json=False):

        params = self._defaults['outputs']

        if to_mongo:
            return self.to_mongo(params['col_reactions'])

        if to_json:
            self.to_json(params['json_reactions'])

    def load_data(self, source='mongo'):
        """
        Overloaded method for loading input data.

        Returns:
            bool: True if successful.
        """

        try:
            params = self._defaults['inputs']

            if source == 'mongo':
                self.side_chains = list(self.from_mongo(params['col_side_chains'], {
                    'type': 'side_chain', 'connection': 'methyl'}))
                self.templates = list(self.from_mongo(params['col_templates'], {'type': 'template'}))
            elif source == 'json':
                self.side_chains = self.from_json(params['json_side_chains'])
                self.templates = self.from_json(params['json_templates'])
            else:
                self.logger.warning(f'The source type \'{source}\' for input data is unrecognized')
        except Exception:
            self.logger.exception('Unexpected exception occured.')
        else:
            return True

        return False

    def generate(self, reactions=[FriedelCafts, TsujiTrost, PictetSpangler, PyrroloIndolene]):

        try:
            for template in self.templates:
                rxn = TemplateOnlyPictetSpangler(Chem.Mol(template['binary']))
                if rxn:
                    data = rxn.data
                    reaction = {'_id': template['_id'] + '_reaction',
                                'type': data.type,
                                'binary': data.binary,
                                'smarts': str(rxn),
                                'rxn_atom_idx': data.rxn_atom_idx,
                                'template': template['_id'],
                                }
                    self.result_data.append(reaction)

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

    def generate_serial(self, reactions=[FriedelCafts, TsujiTrost, PictetSpangler, PyrroloIndolene]):

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

    def generate_from_ids(self, side_chain_ids, template_ids, reactions=[FriedelCafts, TsujiTrost, PictetSpangler, PyrroloIndolene]):

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
            doc = {'_id': side_chain['_id'] + str(chunk + i) + rxn_type[:2],
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
                if atom.GetAtomMapNum() != 0:
                    atom.SetAtomMapNum(0)
                else:
                    atom.SetAtomMapNum(config.SC_WILDCARD_MAP_NUM)
                    utils.atom_to_wildcard(atom)

    def classify_reactions(self, reactions):
        params = self._defaults['template_dependent_rxns']
        dependent_rxns, independent_rxns = [], []
        for func in reactions:
            dependent_rxns.append(func) if func.__name__ in params else independent_rxns.append(func)

        return dependent_rxns, independent_rxns
