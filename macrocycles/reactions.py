
import macrocycles.config as config
import macrocycles.utils as utils
from rdkit import Chem
from collections import namedtuple
from pprint import pprint
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import multiprocessing
from itertools import product

from logging import INFO

LOGGER = utils.create_logger(__name__, INFO)


########################################################################################################################
########################################################################################################################
########################################################################################################################

ReactionInfo = namedtuple('ReactionInfo', 'binary rxn_atom_idx type')


class Reaction():

    def __init__(self, side_chain, template, reacting_atom, name, validation):

        self.validation = validation
        self.name = name
        self.side_chain = side_chain
        self.template = template
        self.reacting_atom = reacting_atom
        self.product = utils.Base.merge(side_chain, template, config.SC_EAS_MAP_NUM,
                                        config.TEMP_EAS_MAP_NUM, clear_map_nums=False) if self else None

    def __repr__(self):
        return 'name: {self.name}\nside_chain: {self.side_chain}\ntemplate: {self.template}\nproduct:{self.product}\n'

    def __str__(self):
        if self.product is None:
            return None

        return '(' + Chem.MolToSmiles(self.side_chain) + '.' + Chem.MolToSmiles(self.template) + ')>>' + Chem.MolToSmiles(self.product)

    def __bool__(self):
        return self.validation(self.reacting_atom)

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


class FriedelCafts(Reaction):

    def __init__(self, side_chain, template, reacting_atom):

        super().__init__(side_chain, template, reacting_atom, 'friedel_crafts', FriedelCafts.is_valid)

    @staticmethod
    def is_valid(atom):

        return atom.GetSymbol() == 'C' and atom.GetTotalNumHs() != 0 and atom.GetIsAromatic()


class TsujiTrost(Reaction):

    def __init__(self, side_chain, template, reacting_atom):
        super().__init__(side_chain, template, reacting_atom, 'tsuji_trost', TsujiTrost.is_valid)

    @staticmethod
    def is_valid(atom):
        return atom.GetSymbol() in ['N', 'O', 'S'] and atom.GetTotalNumHs() > 0


class PictetSpangler(Reaction):

    def __init__(self, side_chain, template, reacting_atom):
        pass

    @staticmethod
    def is_valid(atom):
        pass

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
            self.side_chains = self.from_mongo(params['col_side_chains'], {
                                               'type': 'side_chain', 'connection': 'methyl'})
            self.templates = list(self.from_mongo(params['col_templates'], {'type': 'template'}))
        except Exception:
            self.logger.exception('Unexpected exception occured.')
        else:
            return True

        return False

    def generate(self, reactions=[FriedelCafts, TsujiTrost]):

        try:
            args = product(self.side_chains, self.templates, reactions)
            with multiprocessing.Pool() as pool:
                results = pool.starmap_async(ReactionGenerator.create_reaction, args)
                for rxns, side_chain, template in results.get():
                    self.accumulate_data(rxns, side_chain, template)
        except Exception:
            raise
        else:
            return True

        return False

    def generate_from_ids(self, side_chain_ids, template_ids, reactions=[FriedelCafts, TsujiTrost]):

        try:
            params = self._defaults['inputs']
            self.side_chains = self.from_mongo(params['col_side_chains'], {
                                               'type': 'side_chain', '_id': {'$in': side_chain_ids}})
            self.templates = list(self.from_mongo(params['col_templates'], {
                                  'type': 'template', '_id': {'$in': template_ids}}))
            for side_chain in self.side_chains:
                for template in self.templates:
                    for reaction in reactions:
                        rxn, _, _ = ReactionGenerator.create_reaction(side_chain, template, reaction)
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

        chunk = len(reactions) * self.templates.index(template)
        for i, (smarts, (binary, rxn_atom_idx, rxn_type)) in enumerate(reactions.items()):
            doc = {'_id': side_chain['_id'] + str(chunk + i) + rxn_type[0],
                   'type': rxn_type,
                   'binary': binary,
                   'smarts': smarts,
                   'rxn_atom_idx': rxn_atom_idx,
                   'template': template['_id'],
                   'side_chain': {'_id': side_chain['_id'],
                                  'parent_side_chain': side_chain['parent_side_chain']['_id'],
                                  'conn_atom_idx': side_chain['conn_atom_idx']}}
            self.result_data.append(doc)

    @staticmethod
    def create_reaction(side_chain, template, reaction):

        sc_mol = Chem.Mol(side_chain['binary'])
        temp_mol = Chem.Mol(template['binary'])

        reactions = {}
        ReactionGenerator.tag_wildcard_atom(sc_mol)

        for atom in sc_mol.GetAtoms():

            if atom.GetAtomMapNum() == config.SC_WILDCARD_MAP_NUM:
                continue

            atom.SetAtomMapNum(config.SC_EAS_MAP_NUM)
            rxn = reaction(sc_mol, temp_mol, atom)
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
