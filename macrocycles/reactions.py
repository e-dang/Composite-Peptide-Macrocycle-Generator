
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

    def sc_attach_adjacent_N(self, carboxyl_map_num, nitrogen_map_num, clear_map_nums=True):

        # get nucleophile and backbone
        nucleophile = self.reactants[0]
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

        # create monomer with nucleophile adjacent to nitrogen
        ignored_map_nums = [config.SC_EAS_MAP_NUM, nitrogen_map_num, carboxyl_map_num]
        return utils.Base.merge(nucleophile, backbone, ignored_map_nums=ignored_map_nums, clear_map_nums=clear_map_nums)

########################################################################################################################
########################################################################################################################
########################################################################################################################


class FriedelCafts(Reaction):

    def __init__(self, nucleophile, template, reacting_atom):

        super().__init__(self.is_valid, self.generate_product, 'friedel_crafts', reacting_atom, nucleophile, template)

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

    def __init__(self, nucleophile, template, reacting_atom):

        super().__init__(self.is_valid, self.generate_product, 'tsuji_trost', reacting_atom, nucleophile, template)

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

    def __init__(self, nucleophile, template, reacting_atom):

        super().__init__(self.is_valid, self.generate_product, 'pictet_spangler', reacting_atom, nucleophile, template)

    def is_valid(self):
        if self.valid:
            return True

        nucleophile, template = self.reactants

        # check reacting atom is aromatic carbon with hydrogens
        if not self.reacting_atom.GetIsAromatic() \
                or self.reacting_atom.GetSymbol() != 'C' \
                or self.reacting_atom.GetTotalNumHs() == 0:
            return False

        # check that attachment point of nucleophile is two atoms away from reacting atom
        wild_card = [atom for atom in nucleophile.GetAtoms() if atom.GetAtomicNum() == 0][0]
        paths = Chem.FindAllPathsOfLengthN(nucleophile, 3, useBonds=False, rootedAtAtom=self.reacting_atom.GetIdx())
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
        monomers = self.modify_nucleophile(wild_card)
        template = self.modify_template(substruct_match)
        self.create_reactant(monomers, template)

    def modify_nucleophile(self, wild_card):

        # change nucleophile wild_card back to carbon
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

        # create bond between nucleophile EAS carbon and unmasked aldehyde carbon
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

    def __init__(self, nucleophile, template, reacting_atom):

        super().__init__(self.is_valid, self.generate_product, 'pyrrolo_indolene', reacting_atom, nucleophile, template)

    def is_valid(self):
        if self.valid:
            return True

        nucleophile, template = self.reactants

        # nucleophile must be an indole
        match = nucleophile.GetSubstructMatch(Chem.MolFromSmarts('*c1c[nH]c2ccccc12'))
        if not match:
            return False

        # reaction initiated at carbon that is the attachment point to amino acid backbone, which contains no hydrogens
        if self.reacting_atom.GetTotalNumHs() != 0 \
                or self.reacting_atom.GetSymbol() != 'C' \
                or not self.reacting_atom.GetIsAromatic():
            return False

        # reaction also involves carbon adjacent to reacting_atom which must have a hydrogen
        if 1 not in [atom.GetTotalNumHs() for atom in self.reacting_atom.GetNeighbors()]:
            return False

        # check if carbon is the correct bridged carbon (adjacent to nitrogen)
        wild_card = [atom for atom in nucleophile.GetAtoms() if atom.GetAtomicNum() == 0][0]
        if self.reacting_atom.GetIdx() not in [atom.GetIdx() for atom in wild_card.GetNeighbors()]:
            return False

        self.preprocess()
        self.valid = True
        return True

    def preprocess(self):
        nucleophile = self.modify_nucleophile()
        self.reactants = [nucleophile, Chem.MolFromSmarts(
            f'[*:{config.TEMP_WILDCARD_MAP_NUM}]/C=C/[CH3:{config.TEMP_EAS_MAP_NUM}]')]

    def modify_nucleophile(self):

        # change nucleophile wild card back to carbon
        for atom in self.reactants[0].GetAtoms():
            if atom.GetAtomicNum() == 0:
                atom.SetAtomicNum(6)
                break

        # get nucleophile and backbone
        nucleophile = super().sc_attach_adjacent_N(self.C_TERM_WILDCARD_MAP_NUM, self.NITROGEN_MAP_NUM)

        # add wildcard atom to backbone nitrogen
        unique = set()
        rxn = AllChem.ReactionFromSmarts('[NH2:1][*:2]>>*[NH1:1][*:2]')
        for prod in chain.from_iterable(rxn.RunReactants((nucleophile,))):
            Chem.SanitizeMol(prod)
            unique.add(prod.ToBinary())

        # tag new wildcard atom, peptide nitrogen, and carbon adjacent to reacting atom that has a hydrogen
        nucleophile = Chem.Mol(unique.pop())
        for atom in nucleophile.GetAtoms():
            if atom.GetAtomMapNum() == 0 and atom.GetAtomicNum() == 0:
                atom.SetAtomMapNum(self.N_TERM_WILDCARD_MAP_NUM)
                neighbor = atom.GetNeighbors()[0]
                neighbor.SetAtomMapNum(self.NITROGEN_MAP_NUM)
            elif atom.GetAtomMapNum() == config.SC_EAS_MAP_NUM:
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetTotalNumHs() == 1 and neighbor.GetIsAromatic():
                        neighbor.SetAtomMapNum(self.ADJ_CARBON_MAP_NUM)
                        break

        return nucleophile

    def generate_product(self):

        nucleophile, template = deepcopy(self.reactants)

        atoms = []
        for atom in nucleophile.GetAtoms():
            if atom.GetAtomMapNum() in (config.SC_EAS_MAP_NUM, self.ADJ_CARBON_MAP_NUM):
                atoms.append(atom)
                if len(atoms) == 2:
                    break

        # bond between reacting_atom and adjacent atom becomes single bond
        bond = nucleophile.GetBondBetweenAtoms(atoms[0].GetIdx(), atoms[1].GetIdx())
        bond.SetBondType(Chem.BondType.SINGLE)
        atoms[0].SetNumExplicitHs(atoms[0].GetTotalNumHs() + 1)
        atoms[1].SetNumExplicitHs(atoms[1].GetTotalNumHs() + 1)

        # merge peptide nitrogen to adj carbon
        ignored_map_nums = [self.C_TERM_WILDCARD_MAP_NUM, config.SC_EAS_MAP_NUM, self.N_TERM_WILDCARD_MAP_NUM]
        nucleophile = utils.Base.merge(nucleophile, ignored_map_nums=ignored_map_nums, clear_map_nums=False)

        # merge template with nucleophile
        ignored_map_nums = [self.C_TERM_WILDCARD_MAP_NUM, self.N_TERM_WILDCARD_MAP_NUM,
                            config.TEMP_WILDCARD_MAP_NUM, self.ADJ_CARBON_MAP_NUM, self.NITROGEN_MAP_NUM]
        return utils.Base.merge(nucleophile, template, ignored_map_nums=ignored_map_nums, clear_map_nums=False)

########################################################################################################################
########################################################################################################################
########################################################################################################################


NewReactions = namedtuple('NewReactions', 'reactions nucleophile template')


class ReactionGenerator(utils.Base):
    """
    Class for generating atom mapped reaction SMARTS strings from atom mapped template and side chain SMILES strings.
    Inherits from Base.

    Attributes:
        nucleophiles (list): Contains the different atom mapped side chains.
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
        self.nucleophiles = []
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
            self.nucleophiles = list(self.from_mongo(params['col_nucleophiles'], {
                'type': 'side_chain', 'connection': 'methyl'}))
            self.nucleophiles.extend(list(self.from_mongo(params['col_nucleophiles'], {
                'type': 'monomer', 'group': 'modified_prolines'})))
            self.templates = list(self.from_mongo(params['col_templates'], {'type': 'template'}))
        except Exception:
            self.logger.exception('Unexpected exception occured.')
        else:
            return True

        return False

    def generate(self, reactions=[FriedelCafts, TsujiTrost, PictetSpangler, PyrroloIndolene]):

        try:
            dependent_rxns, independent_rxns = self.classify_reactions(reactions)
            dependent_args = product(self.nucleophiles, dependent_rxns, self.templates)
            independent_args = product(self.nucleophiles, independent_rxns)

            with multiprocessing.Pool() as pool:
                results = pool.starmap_async(ReactionGenerator.create_reaction, dependent_args)
                for rxns, nucleophile, template in results.get():
                    self.accumulate_data(rxns, nucleophile, template)

                results = pool.starmap_async(ReactionGenerator.create_reaction, independent_args)
                for rxns, nucleophile, template in results.get():
                    self.accumulate_data(rxns, nucleophile, template)
        except Exception:
            raise
        else:
            return True

        return False

    def generate_serial(self, reactions=[FriedelCafts, TsujiTrost, PictetSpangler, PyrroloIndolene]):

        try:
            dependent_rxns, independent_rxns = self.classify_reactions(reactions)

            for nucleophile in self.nucleophiles:
                for template in self.templates:
                    for reaction in dependent_rxns:
                        rxn, _, _ = ReactionGenerator.create_reaction(nucleophile, reaction, template)
                        self.accumulate_data(rxn, nucleophile, template)
                for reaction in independent_rxns:
                    rxn, _, template = ReactionGenerator.create_reaction(nucleophile, reaction)
                    self.accumulate_data(rxn, nucleophile, template)
        except Exception:
            raise
        else:
            return True

        return False

    def generate_from_ids(self, nucleophile_ids, template_ids, reactions=[FriedelCafts, TsujiTrost, PictetSpangler, PyrroloIndolene]):

        try:
            params = self._defaults['inputs']
            self.nucleophiles = list(self.from_mongo(params['col_nucleophiles'], {
                'type': 'nucleophile', '_id': {'$in': nucleophile_ids}}))
            self.templates = list(self.from_mongo(params['col_templates'], {
                'type': 'template', '_id': {'$in': template_ids}}))
            dependent_rxns, independent_rxns = self.classify_reactions(reactions)

            for nucleophile in self.nucleophiles:
                for template in self.templates:
                    for reaction in dependent_rxns:
                        rxn, _, _ = ReactionGenerator.create_reaction(nucleophile, reaction, template)
                        self.accumulate_data(rxn, nucleophile, template)

                for reaction in independent_rxns:
                    rxn, _, template = ReactionGenerator.create_reaction(nucleophile, reaction)
                    self.accumulate_data(rxn, nucleophile, template)
        except Exception:
            raise
        else:
            return True

        return False

    def accumulate_data(self, reactions, nucleophile, template):
        """
        Stores all data associated with the different reactions into a dictionary and appends it so self.result_data.

        Args:
            template (dict): The associated data of the templates.
            nucleophile (dict): The associated data of the side chain that the regioisomers are derived from.
            reaction (str): The atom mapped reaction SMARTS string.
        """

        if template is not None:
            chunk = len(reactions) * self.templates.index(template)
            applicable_template = template['_id']
        else:
            chunk = len(reactions)
            applicable_template = 'all'

        for i, (smarts, (binary, rxn_atom_idx, rxn_type)) in enumerate(reactions.items()):
            doc = {'_id': nucleophile['_id'] + str(chunk + i) + rxn_type[:2],
                   'type': rxn_type,
                   'binary': binary,
                   'smarts': smarts,
                   'rxn_atom_idx': rxn_atom_idx,
                   'template': applicable_template,
                   'nucleophile': nucleophile['_id']}
            self.result_data.append(doc)

    @staticmethod
    def create_reaction(nucleophile, reaction, template=None):

        try:
            nuc_mol = Chem.Mol(nucleophile['binary'])
            temp_mol = Chem.Mol(template['binary'])
        except TypeError:
            temp_mol = None

        ReactionGenerator.tag_wildcard_atom(nuc_mol, nucleophile['type'])

        reactions = {}
        for atom in nuc_mol.GetAtoms():

            if atom.GetAtomMapNum() == config.SC_WILDCARD_MAP_NUM:
                continue

            atom.SetAtomMapNum(config.SC_EAS_MAP_NUM)
            rxn = reaction(deepcopy(nuc_mol), deepcopy(temp_mol), atom)
            if rxn:
                reactions[str(rxn)] = rxn.data

            atom.SetAtomMapNum(0)

        return NewReactions(reactions, nucleophile, template)

    @staticmethod
    def tag_wildcard_atom(nucleophile, mol_type):
        if mol_type == 'side_chain':
            matches = nucleophile.GetSubstructMatches(Chem.MolFromSmarts('[CH3]'))
            for pair in matches:
                for atom_idx in pair:
                    atom = nucleophile.GetAtomWithIdx(atom_idx)
                    if atom.GetAtomMapNum() != 0:
                        atom.SetAtomMapNum(0)
                    else:
                        atom.SetAtomMapNum(config.SC_WILDCARD_MAP_NUM)
                        utils.atom_to_wildcard(atom)
        else:
            # change cterm to have wildcard atom
            cterm_rxn = AllChem.ReactionFromSmarts('[*:2][C:1](=O)OH>>[*:2][*:1]')
            mols = {}
            for prod in chain.from_iterable(cterm_rxn.RunReactants((nucleophile,))):
                Chem.SanitizeMol(prod)
                mols[Chem.MolToSmiles(prod)] = prod

            # set atom map num of new wildcard
            mol = list(mols.values())[0]
            for atom in mol.GetAtoms():
                if atom.GetAtomicNum() == 0:
                    atom.SetAtomMapNum(100)

            # change nterm to have wildcard atom
            new_atom = Chem.MolFromSmarts('C')
            new_atom.GetAtomWithIdx(0).SetAtomMapNum(1)
            for atom in mol.GetSubstructMatches(Chem.MolFromSmarts('[NH2,NH1]')):
                atom.SetAtomMapNum(2)
                break

            mol = utils.Base.merge(mol, new_atom, ignored_map_nums=[100])
            for atom in mol.GetAtoms()
                if atom.GetAtomicNum() == 0 and atom.GetAtomicNum() == 0:
                    atom.SetAtomMapNum(101)
                    break

        return atom

    def classify_reactions(self, reactions):
        params = self._defaults['template_dependent_rxns']
        dependent_rxns, independent_rxns = [], []
        for func in reactions:
            dependent_rxns.append(func) if func.__name__ in params else independent_rxns.append(func)

        return dependent_rxns, independent_rxns
