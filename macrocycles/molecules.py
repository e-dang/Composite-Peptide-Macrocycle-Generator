"""
Written by Eric Dang.
github: https://github.com/e-dang
email: edang830@gmail.com
"""

from logging import INFO
from collections import ChainMap, namedtuple, deque, Counter
from copy import deepcopy, copy
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import json
from bson import json_util
from itertools import cycle, chain, islice, product, combinations
from functools import partial
import multiprocessing
from random import sample, seed
from pprint import pprint
import numpy as np
from time import time
import os
import sys

from macrocycles.exceptions import MissingMapNumberError, AtomSearchError, DataNotFoundError, MergeError
from macrocycles.utils import Base, create_logger, read_mols, write_mol, get_user_approval, get_user_atom_idx, atom_to_wildcard, ranges, window
import macrocycles.config as config

LOGGER = create_logger(name=__name__, level=INFO)

########################################################################################################################
########################################################################################################################
########################################################################################################################


class DataInitializer(Base):
    """
    A class for converting data stored in chemdraw files (.sdf) into the proper formats to be operated on. The data is
    then stored in the specified file or database. The data stored in chemdraw files should be the starting point for
    generating macrocyles, such as the templates, parent side chains, backbones, and monomers that won't be derived from
    combinations of side chains and backbones such as natural amino acids and modified prolines. Inherits from Base.
    """

    _defaults = config.DEFAULTS['DataInitializer']

    def __init__(self, logger=LOGGER, make_db_connection=True):
        """
        Constructor.

        Args:
            logger (Logger, optional): The logger. Defaults to LOGGER.
        """

        super().__init__(logger, make_db_connection)

    def initialize_data(self):

        self.load_parent_side_chains()
        self.result_data = []
        self.load_connections()
        self.result_data = []
        self.load_backbones()
        self.result_data = []
        self.load_templates()
        self.result_data = []
        # self.load_rtemplates()
        # self.result_data = []
        self.load_monomers()
        self.result_data = []

    def from_sdf(self, fp, template_doc, func=None):
        """
        Helper function that reads in all molecules in the sdf file defined by self.fp_in, stores the mol's SMILES
        string and kekule SMILES string into the template document, and stores the result in self.result_data.

        Args:
            template_doc (OrderedDict): The dictionary that will contain the SMILES strings.

        Returns:
            bool: True if successful
        """

        try:
            for mol in read_mols(fp):
                doc = deepcopy(template_doc)

                if func is not None:
                    mol, doc = func(mol, doc)

                doc['binary'] = mol.ToBinary()
                Chem.Kekulize(mol)
                doc['kekule'] = Chem.MolToSmiles(mol, kekuleSmiles=True)
                self.result_data.append(doc)
        except (OSError, ValueError, Exception):
            self.logger.exception(f'Variables: fp = {fp}, mol = {Chem.MolToSmiles(mol)}')
            return False

        return True

    def load_parent_side_chains(self, fp=None, group=None, collection=None):
        """
        Top level function that sets up a template document for parent side chains to be stored in, and passes it to
        self.load_files() for filling. If call to self.load_files() is successful then makes call to self.save_data().

        Args:
            group (str, optional): The group name assigned to each molecule in file. Defaults to the last word in the
                file name when split by '_' character.

        Returns:
            bool: True if successful.
        """

        # use defaults or user defined values
        if fp is None:
            fp = self._defaults['inputs']['fp_psc']

        if group is None:
            group = input(f'Enter group name for parent side chains in {fp}: ')

        if collection is None:
            collection = self._defaults['outputs']['col_psc']

        template_doc = {'_id': None,
                        'type': 'parent_side_chain',
                        'binary': None,
                        'kekule': None,
                        'group': group}

        if self.from_sdf(fp, template_doc):
            return self.to_mongo(collection, create_id=True)

        return False

    def load_connections(self, fp=None, collection=None):

        # use defaults or user defined values
        if fp is None:
            fp = self._defaults['inputs']['fp_connections']

        if collection is None:
            collection = self._defaults['outputs']['col_connections']

        template_doc = {'_id': None,
                        'type': 'connection',
                        'binary': None,
                        'kekule': None}
        if self.from_sdf(fp, template_doc, DataInitializer.set_connection_properties):
            return self.to_mongo(collection)

        return False

    @staticmethod
    def set_connection_properties(connection, doc):

        # get _id from user
        Draw.ShowMol(connection)
        doc['_id'] = input('Enter the name of the connection: ')

        # set atom map number of attachment point between parent side chains and connections
        question = 'Enter the index of the atom that will form a bond to the parent side chains: '
        atom_idx = get_user_atom_idx(connection, question)
        connection.GetAtomWithIdx(atom_idx).SetAtomMapNum(config.CONN_MAP_NUM)

        return connection, doc

    def load_backbones(self, fp=None, collection=None):
        """
        Top level function that sets up template documents for the backbone molecules to be stored in and passes it to
        self.load_files() for filling. If call to self.load_files() is successful then makes call to self.save_data().
        Passed in identifier list must be in same order as molecules are read from the sdf file. Default identifiers are
        as follows:
            - 'a' = 'alpha amino acid'
            - 'b' = 'beta2 amino acid'
            - 'c' = 'beta3 amino acid'

        Args:
            identifier (list, optional): A list of strings to be used as IDs. Defaults to ['a', 'b', 'c'].

        Returns:
            bool: True if successful.
        """

        if fp is None:
            fp = self._defaults['inputs']['fp_backbones']

        if collection is None:
            collection = self._defaults['outputs']['col_backbones']

        template_doc = {'_id': None,
                        'type': 'backbone',
                        'binary': None,
                        'kekule': None}

        if self.from_sdf(fp, template_doc, DataInitializer.set_backbone_properties):
            return self.to_mongo(collection)

        return False

    @staticmethod
    def set_backbone_properties(backbone, doc):

        # get _id from user
        Draw.ShowMol(backbone)
        doc['_id'] = input('Enter the name of the backbone: ')

        # get atom index from user
        question = 'Enter the index of the atom that will form a bond with the side chains: '
        atom_idx = get_user_atom_idx(backbone, question)
        backbone.GetAtomWithIdx(atom_idx).SetAtomMapNum(config.BB_MAP_NUM)

        return backbone, doc

    def load_monomers(self, fp=None, group=None, backbone=None, required=None, collection=None):
        """
        Top level function that sets up template document for monomers (that are not derived from modified side chains
        and backbones) to be stored in and passes it to self.load_files() for filling. If call to self.load_files() is
        successful then makes call to self.save_data().

        Args:
            group (str, optional): The group name assigned to each molecule in file. Defaults to the last word in the
                file name when split by '_' character.
            required (bool, optional): The value of 'required' stored in the documents. This value determines if the
                peptides these monomers are present in are valid if no other required monomers are present. Defaults to
                False.
            backbone (str, optional): The backbone type these monomers are made of. Defaults to None.

        Returns:
            bool: True if successful.
        """

        # use defaults or user defined values
        if fp is None:
            fp = self._defaults['inputs']['fp_monomers']

        if backbone is None:
            backbone = input('Enter the backbone type these monomers are made of: ')
            bb_doc = self.from_mongo(self._defaults['outputs']['col_backbones'], {'_id': backbone})[0]
            if bb_doc is None:
                raise DataNotFoundError(f'The backbone with _id \'{backbone}\' could not be found.')

        if group is None:
            group = input(f'Enter group name for monomers in file {fp}: ')

        if collection is None:
            collection = self._defaults['outputs']['col_monomers']

        template_doc = {'_id': None,
                        'type': 'monomer',
                        'binary': None,
                        'kekule': None,
                        'backbone': {'_id': bb_doc['_id'],
                                     'binary': bb_doc['binary']},
                        'side_chain': None,
                        'group': group,
                        'required': required}

        if required is None and self.from_sdf(fp, template_doc, DataInitializer.set_monomer_properties):
            return self.to_mongo(collection, create_id=True)

        if self.from_sdf(fp, template_doc):
            return self.to_mongo(collection, create_id=True)

        return False

    @staticmethod
    def set_monomer_properties(monomer, doc):

        # get required field from user
        Draw.ShowMol(monomer)
        doc['required'] = get_user_approval('Should this monomer be in the set of required monomers (y/n)?: ')

        return monomer, doc

    def load_templates(self, fp=None, collection=None):
        """
        Top level function that sets up template documents for the template molecules to be stored in and passes it to
        self.load_files() for filling. If call to self.load_files() is successful then makes call to self.save_data().
        Passed in identifier list must be in same order as molecules are read from the sdf file. Defaults indentifiers
        are as follows:
            - 't1' = template1
            - 't2' = template2
            - 't3' = template3

        Args:
            ID (list, optional): A list of strings to be used as IDs. Defaults to ['t1', 't2', 't3'].

        Returns:
            bool: True if successful.
        """

        # use defaults or user defined values
        if fp is None:
            fp = self._defaults['inputs']['fp_templates']

        if collection is None:
            col_temp = self._defaults['outputs']['col_templates']
            col_rtemp = self._defaults['outputs']['col_rtemplates']

        temp_doc = {'_id': None,
                        'type': 'template',
                        'binary': None,
                        'kekule': None,
                        'original': None}

        rtemp_doc = {'_id': None,
                        'type': 'template',
                        'binary': None,
                        'smarts': None,
                        'original': None}

        temp_data, rtemp_data = [], []
        for mol in read_mols(fp):
            doc = deepcopy(temp_doc)
            rdoc = deepcopy(rtemp_doc)

            mol, doc = DataInitializer.set_template_properties(mol, doc)
            rdoc['_id'] = doc['_id']
            rdoc['original'] = doc['original']
            rtemp_data.append(DataInitializer.modify_rtemplate(deepcopy(mol), rdoc))
            temp_data.append(DataInitializer.modify_template(deepcopy(mol), doc))

        self.result_data = temp_data
        temp_result = self.to_mongo(col_temp)
        self.result_data = rtemp_data
        rtemp_result = self.to_mongo(col_rtemp)
        return temp_result and rtemp_result

    @staticmethod
    def set_template_properties(template, doc):

        # get SMILES of original template before modification
        doc['original'] = Chem.MolToSmiles(template)

        # get _id from user
        Draw.ShowMol(template)
        doc['_id'] = input('Enter the name of the template: ')

        # set atom map number for connecting atom between template and peptide
        question = 'Enter the index of the carbon atom that will form amide linkage to the peptides: '
        substruct = DataInitializer.generalize_substruct(template, config.TEMP_WILDCARD_MAP_NUM, question)
        template = Chem.DeleteSubstructs(template, Chem.MolFromSmiles(substruct))

        # remove leaving group at EAS reaction site
        question = 'Enter the index of the carbon atom that will close the macrocycle ring through an EAS reaction: '
        substruct = DataInitializer.generalize_substruct(template, config.TEMP_EAS_MAP_NUM, question)
        template = Chem.DeleteSubstructs(template, Chem.MolFromSmiles(substruct))

        return template, doc

    @staticmethod
    def modify_rtemplate(rtemplate, doc):

        mol = Chem.RWMol(rtemplate)
        for atom in mol.GetAtoms():
            if atom.GetAtomMapNum() == config.TEMP_WILDCARD_MAP_NUM:
                atom_to_wildcard(atom)
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() == 'O':
                        break
        mol.RemoveAtom(neighbor.GetIdx())
        Chem.SanitizeMol(mol)
        mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
        doc['binary'] = mol.ToBinary()
        Chem.Kekulize(mol)
        doc['smarts'] = Chem.MolToSmiles(mol, kekuleSmiles=True)

        return doc

    @staticmethod
    def modify_template(template, doc):

        for atom in template.GetAtoms():
            if atom.GetAtomMapNum() == config.TEMP_EAS_MAP_NUM:
                atom.SetAtomMapNum(0)
        doc['binary'] = template.ToBinary()
        Chem.Kekulize(template)
        doc['kekule'] = Chem.MolToSmiles(template, kekuleSmiles=True)

        return doc

    def load_rtemplates(self, fp=None, collection=None):

        if fp is None:
            fp = self._defaults['inputs']['fp_templates']

        if collection is None:
            collection = self._defaults['outputs']['col_rtemplates']

        template_doc = {'_id': None,
                        'type': 'template',
                        'binary': None,
                        'smarts': None,
                        'original': None}

        if self.from_sdf(fp, template_doc, DataInitializer.set_rtemplate_properties):
            for doc in self.result_data:
                doc['smarts'] = doc['kekule']
                del doc['kekule']
            return self.to_mongo(collection)

        return False

    @staticmethod
    def set_rtemplate_properties(template, doc):

        # get SMILES of original template before modification
        doc['original'] = Chem.MolToSmiles(template)

        # get _id from user
        Draw.ShowMol(template)
        doc['_id'] = input('Enter the name of the template: ')

        # set atom map number of the EAS reacting atom and get leaving group
        question = 'Enter the index of the atom that will under go an EAS reaction to close the macrocycle ring: '
        substruct = DataInitializer.generalize_substruct(template, config.TEMP_EAS_MAP_NUM, question)
        template = Chem.DeleteSubstructs(template, Chem.MolFromSmiles(substruct))

        # change desired atom to wildcard and set atom map number
        question = 'Enter index of the atom that will be turned into a wildcard atom: '
        substruct = DataInitializer.generalize_substruct(template, config.TEMP_WILDCARD_MAP_NUM, question)
        template = Chem.DeleteSubstructs(template, Chem.MolFromSmiles(substruct))
        for atom in template.GetAtoms():
            if atom.GetAtomMapNum() == config.TEMP_WILDCARD_MAP_NUM:
                atom_to_wildcard(atom)

        return template, doc

    @staticmethod
    def generalize_substruct(mol, map_num, question):

        while True:
            mol_copy = deepcopy(mol)
            atom_idx = get_user_atom_idx(mol_copy, question)
            atom = mol_copy.GetAtomWithIdx(atom_idx)
            atom.SetAtomMapNum(map_num)

            for bond in atom.GetBonds():
                fragmented_mol = Chem.FragmentOnBonds(mol_copy, [bond.GetIdx()], dummyLabels=[(0, 0)])
                for frag in Chem.GetMolFrags(fragmented_mol, asMols=True):
                    if all([atom.GetAtomMapNum() == 0 for atom in frag.GetAtoms()]):
                        Draw.ShowMol(frag)
                        if get_user_approval('Is this the substructure to generalize (y/n)?: '):
                            mol.GetAtomWithIdx(atom_idx).SetAtomMapNum(map_num)
                            return Chem.MolToSmiles(frag, isomericSmiles=False).strip('*')

            print('Please select a substructure')

    def to_json(self):

        params = self._defaults['outputs']

        parent_side_chains = list(self.from_mongo(params['col_psc'], {'type': 'parent_side_chain'}))
        connections = list(self.from_mongo(params['col_connections'], {'type': 'connection'}))
        backbones = list(self.from_mongo(params['col_backbones'], {'type': 'backbone'}))
        monomers = list(self.from_mongo(params['col_monomers'], {
                        'type': 'monomer', 'side_chain': {'$not': {'$type': 'object'}}}))
        templates = list(self.from_mongo(params['col_templates'], {'type': 'template'}))
        rtemplates = list(self.from_mongo(params['col_rtemplates'], {'type': 'template'}))

        # TODO: Catch errrors properly
        with open(params['fp_psc'], 'w') as file:
            json.dump(json.loads(json_util.dumps(parent_side_chains)), file)

        with open(params['fp_connections'], 'w') as file:
            json.dump(json.loads(json_util.dumps(connections)), file)

        with open(params['fp_backbones'], 'w') as file:
            json.dump(json.loads(json_util.dumps(backbones)), file)

        with open(params['fp_monomers'], 'w') as file:
            json.dump(json.loads(json_util.dumps(monomers)), file)

        with open(params['fp_templates'], 'w') as file:
            json.dump(json.loads(json_util.dumps(templates)), file)

        with open(params['fp_rtemplates'], 'w') as file:
            json.dump(json.loads(json_util.dumps(rtemplates)), file)

    def from_json(self):

        params = self._defaults['outputs']
        try:
            with open(params['fp_psc'], 'r') as file:
                self.result_data = json_util.loads(json_util.dumps(json.load(file)))
                self.to_mongo(config.COL1)

            with open(params['fp_connections'], 'r') as file:
                self.result_data = json_util.loads(json_util.dumps(json.load(file)))
                self.to_mongo(config.COL1)

            with open(params['fp_backbones'], 'r') as file:
                self.result_data = json_util.loads(json_util.dumps(json.load(file)))
                self.to_mongo(config.COL1)

            with open(params['fp_monomers'], 'r') as file:
                self.result_data = json_util.loads(json_util.dumps(json.load(file)))
                self.to_mongo(config.COL1)

            with open(params['fp_templates'], 'r') as file:
                self.result_data = json_util.loads(json_util.dumps(json.load(file)))
                self.to_mongo(config.COL1)

            with open(params['fp_rtemplates'], 'r') as file:
                self.result_data = json_util.loads(json_util.dumps(json.load(file)))
                self.to_mongo(config.COL2)
        except (OSError, json.JSONDecodeError):  # TODO: Catch errors properly
            self.logger.exception(f'Failed to load the json file(s)')
            return False

        return True

########################################################################################################################
########################################################################################################################
########################################################################################################################

SideChainInfo = namedtuple('SideChainInfo', 'smiles conn_atom_idx')
NewSideChains = namedtuple('NewSideChains', 'side_chain parent_sc connection')


class SideChainGenerator(Base):
    """
    Class for attaching varying length alkyl chains / attachment points to side chain. Inherits from Base.

    Attributes:
        parent_side_chains (list): Contains the parent_side_chains and associated data as dictionaries.
        connections (list): Contains the atom mapped SMARTS strings of the alkyl attachment chain and the corresponding
            modification array as a Connections named tuple.
    """

    _defaults = config.DEFAULTS['SideChainGenerator']

    def __init__(self, logger=LOGGER, make_db_connection=True):
        """
        Initializer.

        Args:
            logger (Logger)
        """

        # I/O
        super().__init__(logger, make_db_connection)

        # data
        self.parent_side_chains = []
        self.connections = []

    def save_data(self):

        params = self._defaults['outputs']

        return self.to_mongo(params['col_side_chains'])

    def load_data(self, **kwargs):
        """
        Overloaded method for loading input data using default locations defined in config.py. Tries to load from the
        database first, but if connection is not established then tries to load from a json file.

        Returns:
            bool: True if successful.
        """

        params = ChainMap(kwargs, self._defaults['inputs'])

        try:
            self.parent_side_chains = self.from_mongo(params['col_psc'], {'type': 'parent_side_chain'})
            self.connections = list(self.from_mongo(params['col_connections'], {'type': 'connection'}))
        except Exception:
            self.logger.exception('Unexpected exception occured.')
        else:
            return True

        return False

    def generate(self):
        """
        Main driver of class functionality. Calls self.alternate_connection_point() on each
        molecule in self.parent_side_chains with each connection type in self.connections and calls
        self.accumulate_mols() on the resulting data.

        Returns:
            bool: True if successful.
        """

        try:
            args = product(self.parent_side_chains, self.connections)
            with multiprocessing.Pool() as pool:
                results = pool.starmap_async(SideChainGenerator.create_side_chain, args)
                for side_chain, parent_sc, connection in results.get():
                    self.accumulate_data(side_chain, parent_sc, connection)
        except MissingMapNumberError:
            self.logger.exception(f'connection: {connection}')
        except Exception:
            raise
        else:
            return True

        return False

    def generate_from_ids(self, parent_sc_ids, connection_ids, **kwargs):

        params = ChainMap(kwargs, self._defaults['inputs'])

        try:
            self.parent_side_chains = list(self.from_mongo(params['col_psc'], {'_id': {'$in': parent_sc_ids}}))
            self.connections = list(self.from_mongo(params['col_connections'], {'_id': {'$in': connection_ids}}))

            results = []
            for parent_sc in self.parent_side_chains:
                for connection in self.connections:
                    side_chains, _, _ = SideChainGenerator.create_side_chain(parent_sc, connection, self.logger)
                    results.append(side_chains)
                    self.accumulate_data(side_chains, parent_sc, connection)
        except MissingMapNumberError:
            self.logger.exception(f'connection: {connection}')
        else:
            self.logger.info(f'Successfully generated {len(self.result_data)} side chains')
            return results

        return False

    def accumulate_data(self, side_chains, parent_sc, connection):
        """
        Stores all data associated with the modified side chain in a dictionary and appends it to self.result_data.

        Args:
            unique_mols (iterable): Contains the unique SMILES strings of the modified side chain.
            parent (dict): The assocaited data of parent side chain from which the modified side chain was derived.
            modifications (list): A list containing ids of the modifications that were made to the parent side chain.
        """

        chunk = len(side_chains) * self.connections.index(connection)
        for i, (binary, (kekule, conn_atom_idx)) in enumerate(side_chains.items()):
            doc = {'_id': parent_sc['_id'] + str(chunk + i),
                   'type': 'side_chain',
                   'binary': binary,
                   'kekule': kekule,
                   'conn_atom_idx': conn_atom_idx,
                   'connection': connection['_id'],
                   'parent_side_chain': {'_id': parent_sc['_id'],
                                        'group': parent_sc['group']}}
            self.result_data.append(doc)

    @staticmethod
    def create_side_chain(parent_sc, connection, logger=None):
        """
        Creates a set of new molecules by attaching an alkyl chain (which becomes the attachment point to peptide
        backbone) to every eligble position on the side chain. Eligiblity of an atom is defined as:
            Carbon - Must have 1 or 2 hydrogens
            Nitrogen, Oxygen, Sulfur - Must have > 0 hydrogens

        Args:
            mol (str): The SMILES string of the parent side chain molecule.
            connection_tup (Connections): A namedtuple containing the atom mapped alkyl attachment chain and
                modifications array.

        Raises:
            MissingMapNumberError: If connection is not atom mapped.

        Returns:
            dict: Containing the unique side chain SMILES strings as keys and the corresponding kekule SMILES and
                modifications as values.
        """

        conn_mol = Chem.Mol(connection['binary'])
        parent_mol = Chem.Mol(parent_sc['binary'])

        # check if connecting atom is atom mapped
        SideChainGenerator.is_valid_connection(conn_mol)

        # make attachment at each atom
        unique_mols = {}
        for atom in parent_mol.GetAtoms():

            # detetmine atom eligibility
            if SideChainGenerator.is_valid_atom(atom):
                atom.SetAtomMapNum(config.PSC_MAP_NUM)
            else:
                continue

            # merge parent side chain with conenction and record results
            try:
                # side_chain = Base.merge(parent_mol, conn_mol, config.PSC_MAP_NUM, config.CONN_MAP_NUM)
                side_chain = Base.merge(parent_mol, conn_mol)
            except (ValueError, MergeError) as error:
                if logger:
                    logger.exception(f'Parent Side Chain: {Chem.MolToSmiles(parent_mol)}')
                else:
                    print(error)
            else:
                atom.SetAtomMapNum(0)
                binary = side_chain.ToBinary()
                Chem.Kekulize(side_chain)
                unique_mols[binary] = SideChainInfo(Chem.MolToSmiles(side_chain, kekuleSmiles=True), atom.GetIdx())

        return NewSideChains(unique_mols, parent_sc, connection)

    @staticmethod
    def is_valid_atom(atom):

        valid_carbon = atom.GetSymbol() == 'C' and 0 < atom.GetTotalNumHs() < 3
        valid_hetero = atom.GetSymbol() in ('N', 'O', 'S') and atom.GetTotalNumHs() != 0

        return valid_carbon or valid_hetero

    @staticmethod
    def is_valid_connection(connection):

        if config.CONN_MAP_NUM not in [atom.GetAtomMapNum() for atom in connection.GetAtoms()]:
            raise MissingMapNumberError('Connection molecule missing atom map number')

        return True

########################################################################################################################
########################################################################################################################
########################################################################################################################


NewMonomer = namedtuple('NewMonomer', 'monomer side_chain backbone stereo')


class MonomerGenerator(Base):
    """
    Class for combining side_chains and amino acid backbones to form monomers. Inherits from Base.

    Attributes:
        backbones (list): Contains the different amino acid backbones.
        side_chains (list): Contains the side_chains and the associated data as dictionaries.
    """

    _defaults = config.DEFAULTS['MonomerGenerator']

    def __init__(self, logger=LOGGER, make_db_connection=True):
        """
        Initializer.

        Args:
            required (bool): The value of the 'required' field to be stored in the monomer documents.
        """

        # I/O
        super().__init__(logger, make_db_connection)

        # data
        self.side_chains = []
        self.backbones = []

    def save_data(self):

        params = self._defaults['outputs']

        return self.to_mongo(params['col_monomers'])

    def load_data(self, **kwargs):
        """
        Overloaded method for loading input data.

        Returns:
            bool: True if successful.
        """

        params = ChainMap(kwargs, self._defaults['inputs'])

        try:
            self.side_chains = self.from_mongo(params['col_side_chains'], {'type': 'side_chain'})
            self.backbones = list(self.from_mongo(params['col_backbones'], {'type': 'backbone'}))
        except Exception:
            self.logger.exception('Unexcepted exception occured.')
        else:
            return True

        return False

    def generate(self, stereochem=['CW', 'CCW'], required=False):

        try:
            args = product(self.side_chains, self.backbones, stereochem)
            with multiprocessing.Pool() as pool:
                results = pool.starmap_async(MonomerGenerator.create_monomer, args)
                for monomer, side_chain, backbone, stereo in results.get():
                    self.accumulate_data(monomer, side_chain, backbone, required, stereo)
        except AtomSearchError:
            self.logger.exception('')
        except Exception:
            raise
        else:
            return True

        return False

    def generate_from_ids(self, side_chain_ids, backbone_ids, stereochem=['CW', 'CCW'], required=False, **kwargs):
        """
        Top level function that calls self.create_monomer() on each side_chain - backbone - stereochemistry combination
        and calls self.accumulate_mols() with the resulting monomers.

        Returns:
            bool: True if successful.
        """

        params = ChainMap(kwargs, self._defaults['inputs'])

        try:
            self.side_chains = list(self.from_mongo(params['col_side_chains'], {'_id': {'$in': side_chain_ids}}))
            self.backbones = list(self.from_mongo(params['col_backbones'], {'_id': {'$in': backbone_ids}}))

            results = []
            for side_chain in self.side_chains:
                for backbone in self.backbones:
                    for stereo in stereochem:
                        monomer, _, _, _ = MonomerGenerator.create_monomer(side_chain, backbone, stereo, self.logger)
                        results.append(monomer)
                        self.accumulate_data(monomer, side_chain, backbone, required, stereo)
        except AtomSearchError:
            self.logger.exception('')
        except Exception:
            self.logger.exception('Unexpected exception occured.')
        else:
            self.logger.info(f'Successfully generated {len(self.result_data)} monomers')
            return results

        return False

    def accumulate_data(self, monomers, side_chain, backbone, required, stereo):
        """
        Stores all data associated with the monomer in a dictionary and appends it to self.result_data.

        Args:
            monomers (dict): The dictionary containing the monomers' SMILES string as keys and kekule SMILES as values.
            backbone (dict): The dictionary containing the associated backbone data.
            side_chain (dict): The dictionary containing the associated side chain data.
        """

        bb_ind = self.backbones.index(backbone)
        stereo_ind = 'f' if stereo == 'CW' else 'r'
        comb_ind = stereo_ind + str(bb_ind)

        for binary, kekule in monomers.items():
            doc = {'_id': comb_ind + side_chain['_id'].upper(),
                   'type': 'monomer',
                   'binary': binary,
                   'kekule': kekule,
                   'required': required,
                   'group': side_chain['parent_side_chain']['group'],
                   'backbone': {'_id': backbone['_id'],
                                'binary': backbone['binary']},
                   'side_chain': {'_id': side_chain['_id'],
                                  'parent_side_chain': side_chain['parent_side_chain']['_id'],
                                  'conn_atom_idx': side_chain['conn_atom_idx']}}
            self.result_data.append(doc)

    @staticmethod
    def create_monomer(side_chain, backbone, stereo, logger=None):
        """
        Connects the side chain to the backbone at the designated attachment points.

        Args:
            backbone (rdkit Mol): The backbone to which the side chain will be attached to.
            side_chain (rdkit Mol): The side chain being attached to the backbone structure.
            stereo (str): The stereochemistry at the attachment point of the backbone and side chain.

        Returns:
            dict: A dictionary containing the monomers' SMILES string as keys and kekule SMILES as values.
        """

        sc_mol = Chem.Mol(side_chain['binary'])
        bb_mol = Chem.Mol(backbone['binary'])

        # find attachment point(s) at terminal end of alkyl chain
        matches = sc_mol.GetSubstructMatches(Chem.MolFromSmarts('[CH3]'))
        if not len(matches):
            raise AtomSearchError(f'Couldn\'t find attchment point of: {Chem.MolToSmiles(sc_mol)}')

        monomers = {}
        for pairs in matches:

            # set atom map number for attachment point of side chain
            for atom_idx in pairs:
                atom = sc_mol.GetAtomWithIdx(atom_idx)
                atom.SetAtomMapNum(config.SC_MAP_NUM)
                break

            # connect monomer and backbone
            try:
                monomer = Base.merge(sc_mol, bb_mol, stereo=stereo)
            except ValueError:
                if logger:
                    logger.exception(f'Sanitize Error! Side Chain: {Chem.MolToSmiles(sc_mol)}, '
                                     f'backbone: {Chem.MolToSmiles(bb_mol)}')
            else:
                atom.SetAtomMapNum(0)
                binary = monomer.ToBinary()
                Chem.Kekulize(monomer)
                monomers[binary] = Chem.MolToSmiles(monomer, kekuleSmiles=True)

        return NewMonomer(monomers, side_chain, backbone, stereo)

########################################################################################################################
########################################################################################################################
########################################################################################################################


NewPeptide = namedtuple('NewPeptide', 'peptide monomers')


class PeptideGenerator(Base):
    """
    Class for joining monomers together to form a peptide. Inherits from Base.

    Attributes:
        monomers (list): Contains the monomers and the associated data.
    """

    _defaults = config.DEFAULTS['PeptideGenerator']

    def __init__(self, logger=LOGGER, make_db_connection=True):
        """
        Initializer.

        Args:
        """

        # I/O
        super().__init__(logger, make_db_connection)

        # data
        self.monomers = []

    def save_data(self):

        params = self._defaults['outputs']

        return self.to_mongo(params['col_peptides'])

    def load_data(self):
        """
        Overloaded method for loading input data.

        Returns:
            bool: True if successful.
        """

        params = self._defaults['inputs']

        try:
            self.monomers = self.from_mongo(params['col_monomers'], {'type': 'monomer'})
        except Exception:
            self.logger.exception('Unexcepted exception occured.')
        else:
            return True

        return False

    def generate(self, length, num_peptides=None, random=False, num_jobs=1, job_id=1):
        """
        Top level function that generates the cartesian product of monomers in self.monomers. create_peptide() is then
        called in a parallel fashion on each tuple of monomers is the cartesian product. The resulting peptide and
        monomers are then passed to self.accumulate_data(). If self.result_data is about to exceed the capacity of a
        python list, an automatic call to self.save_data() is made and self.result_data is cleared.

        Args:
            length (int): The number of monomers per peptide.
            num_peptides (int, optional): The number of peptides to make. If None, then all peptides in the cartesian
                product as produced. Defaults to None.
            num_jobs (int, optional): The number of jobs submitted in a job array. Used for running this method on a
                super computer. Does not do anything if num_peptides is specified. Defaults to 1.
            job_id (int, optional): The job ID of the specific job in the job array running this function. Used for
                running this method on a super computer. Does not do anything if num_peptides is specified. Defaults
                to 1.

        Returns:
            bool: True if successful.
        """

        try:
            # determine how to run method
            if num_peptides:
                if random:
                    self.monomers = list(self.monomers)
                    seed(time())
                    monomers = [sample(self.monomers, length) for i in range(num_peptides)]
                else:
                    monomers = islice(product(self.monomers, repeat=length), num_peptides)
            else:
                start, stop = ranges(len(self.monomers) ** length, num_jobs)[job_id - 1]
                monomers = islice(product(self.monomers, repeat=length), start, stop)

            # perform peptide generation in parallel
            invalids = 0
            with multiprocessing.Pool() as pool:
                result = pool.map_async(PeptideGenerator.create_peptide, monomers)
                for peptide, monomers in result.get():
                    if None in (peptide, monomers):
                        invalids += 1
                    else:
                        self.accumulate_data(peptide, monomers)

                    # save data if reaching size limit of python lists
                    if len(self.result_data) > config.CAPACITY:
                        self.save_data()
                        self.result_data = []
        except Exception:
            self.logger.exception('Unexpected exception occured.')
        else:
            self.logger.info(f'Successfully created {len(self.result_data)} valid peptides and discarded {invalids} '
                            'invalid peptides!')
            return True

        return False

    def generate_from_ids(self, *monomer_ids):

        params = self._defaults['inputs']

        try:
            count = 0
            results = []
            for monomers in monomer_ids:
                counts = Counter(monomers)
                monomers = self.from_mongo(params['col_monomers'], {'type': 'monomer', '_id': {'$in': monomers}})
                monomers = [monomer for monomer in monomers for i in range(counts[monomer['_id']])]
                peptide, monomers = PeptideGenerator.create_peptide(monomers)
                if None in (peptide, monomers):
                    count = 0
                else:
                    results.append(peptide)
                    self.accumulate_data(peptide, monomers)
        except Exception:
            raise
        else:
            return results, count

        return False

    def accumulate_data(self, peptide, monomers):
        """
        Stores all data associated with the peptide in a dictionary and appends it to self.result_data.

        Args:
            peptide (dict): The dictionary containing the peptide's SMILES string as keys and a list containing the
                kekule SMILES and the N-terminal monomer's backbone type as values.
            monomers (dict): The dictionary containing the associated monomer data.
        """

        monomer_data = [{key: value for key, value in monomer.items() if key in ('_id', 'side_chain')} for monomer in monomers]
        pep_id = ''.join([monomer['_id'] for monomer in monomer_data])

        for binary, kekule in peptide.items():
            doc = {'_id': pep_id,
                    'type': 'peptide',
                    'binary': binary,
                    'kekule': kekule,
                    'monomers': monomer_data}
            self.result_data.append(doc)

    @staticmethod
    def create_peptide(monomers, logger=None):
        """
        Static class method. Connects monomers together creating a peptide and verifies that the peptide is composed of
        at least one required monomer.

        Args:
            monomers (list): Contains the monomer documents.
            logger (Logger, optional): A logger used to log errors in the method. Defaults to None.

        Returns:
            dict: A dictionary containing the peptide SMILES string as a key and a list containing the kekule SMILES
                string and N-terminal monomer's backbone type as a value.
            list: The list of monomer documents used to create the peptide.
        """

        # check that peptide is valid before constructing
        if not PeptideGenerator.is_valid_monomers(monomers):
            return None, None

        # begin conneting each monomer in monomers
        for i, monomer in enumerate(monomers):

            monomer_mol = Chem.Mol(monomer['binary'])
            backbone = Chem.Mol(monomer['backbone']['binary'])

            # start peptide with first monomer
            if i == 0:
                peptide = monomer_mol
                backbone_prev = backbone
                continue

            # assign atom map numbers
            monomer_old_attach = PeptideGenerator.tag_monomer_n_term(monomer_mol, backbone)
            carboxyl_atom, pep_old_attach = PeptideGenerator.tag_peptide_c_term(peptide, backbone_prev)

            # remove oxygen atom from carboxyl
            peptide = Chem.RWMol(peptide)
            peptide.RemoveAtom(carboxyl_atom)

            # connect peptide and monomer
            try:
                peptide = Base.merge(peptide, monomer_mol)
            except ValueError:
                if logger:
                    logger.exception(f'Sanitize Error! monomer = {Chem.MolToSmiles(monomer_mol)}, '
                                     f'peptide = {Chem.MolToSmiles(peptide)}')
                else:
                    print(f'Sanitize Error! monomer = {Chem.MolToSmiles(monomer_mol)}, '
                          f'peptide = {Chem.MolToSmiles(peptide)}')
            else:
                pep_old_attach.SetAtomMapNum(0)
                monomer_old_attach.SetAtomMapNum(0)
                backbone_prev = backbone

        # structure return data
        binary = peptide.ToBinary()
        Chem.Kekulize(peptide)
        result = {binary: Chem.MolToSmiles(peptide, kekuleSmiles=True)}
        return NewPeptide(result, monomers)

    @staticmethod
    def is_valid_monomers(monomers):
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

    @staticmethod
    def tag_monomer_n_term(monomer, backbone):

        matches = monomer.GetSubstructMatches(backbone)
        for pair in matches:
            for atom_idx in pair:
                atom = monomer.GetAtomWithIdx(atom_idx)
                if atom.GetSymbol() == 'N' and atom.GetTotalNumHs() != 0:
                    atom.SetAtomMapNum(config.MONO_NITROGEN_MAP_NUM)
                    return atom

    @staticmethod
    def tag_peptide_c_term(peptide, backbone):

        matches = peptide.GetSubstructMatches(backbone)
        for pair in matches:
            for atom_idx in pair:
                atom = peptide.GetAtomWithIdx(atom_idx)
                if atom.GetSymbol() == 'O' and atom.GetTotalNumHs() == 1:
                    carboxyl_atom = atom_idx
                elif atom.GetSymbol() == 'C' and atom.GetTotalNumHs() == 0 and \
                        atom.GetHybridization() == Chem.rdchem.HybridizationType.SP2:
                    atom.SetAtomMapNum(config.PEP_CARBON_MAP_NUM)
                    attachment_atom = atom

        return carboxyl_atom, attachment_atom

########################################################################################################################
########################################################################################################################
########################################################################################################################

NewTPHybrid = namedtuple('NewTPHybrid', 'tp_hybrid peptide template')

class TPHybridGenerator(Base):
    """
    Class for joining templates and peptides together. Inherits from Base.

    Attributes:
        templates (list): Contains the templates and the associated data.
        peptides (list): Contains the peptides and the associated data.
    """

    _defaults = config.DEFAULTS['TPHybridGenerator']

    def __init__(self, logger=LOGGER, make_db_connection=True):
        """
        Initializer.

        Args:
        """

        # I/O
        super().__init__(logger, make_db_connection)

        # data
        self.peptides = []
        self.templates = []

    def save_data(self):

        params = self._defaults['outputs']

        return self.to_mongo(params['col_tp_hybrids'])

    def load_data(self):
        """
        Overloaded method for loading input data.

        Returns:
            bool: True if successful.
        """

        params = self._defaults['inputs']
        try:
            self.peptides = self.from_mongo(params['col_peptides'], {'type': 'peptide'})
            self.templates = list(self.from_mongo(params['col_templates'], {'type': 'template'}))
        except Exception:
            self.logger.exception('Unexcepted exception occured.')
        else:
            return True

        return False

    def generate(self):

        try:
            args = product(self.peptides, self.templates)
            with multiprocessing.Pool() as pool:
                results = pool.starmap_async(TPHybridGenerator.create_tp_hybrid, args)
                for tp_hybrid, peptide, template in results.get():
                    self.accumulate_data(tp_hybrid, peptide, template)
        except Exception:
            raise
        else:
            return True

        return False

    def generate_from_ids(self, peptide_ids, template_ids):
        """
        Main driver function that calls self.modify_templates() on each template, and self.connect_template_peptide() on
        each unique peptide-template pair. The resulting tp_hybrid is then passed to self.accumulate_data().

        Returns:
            bool: True if successful.
        """


        try:
            params = self._defaults['inputs']
            self.peptides = self.from_mongo(params['col_peptides'], {'type': 'peptide', '_id': {'$in': peptide_ids}})
            self.templates = list(self.from_mongo(params['col_templates'],
                                                  {'type': 'template', '_id': {'$in': template_ids}}))

            results = []
            for peptide in self.peptides:
                for template in self.templates:
                    tp_hybrid, _, _ = TPHybridGenerator.create_tp_hybrid(peptide, template, self.logger)
                    self.accumulate_data(tp_hybrid, peptide, template)
                    results.append(tp_hybrid)
        except Exception:
            raise
        else:
            return results

        return False

    def accumulate_data(self, tp_hybrid, peptide, template):
        """
        Stores all data associated with the tp_hybrid in a dictionary and appends it to self.result_data.

        Args:
            tp_hybrid (dict): The dictionary containing the tp_hybrid SMILES string as keys and kekule SMILES as values.
            template (dict): The dictionary containing the associated template data.
            peptide (dict): The dictionary containing the associated peptide data.
        """

        chunk = len(tp_hybrid) * self.templates.index(template)

        for i, (binary, kekule) in enumerate(tp_hybrid.items()):
            doc = {'_id': peptide['_id'] + str(chunk + i),
                    'type': 'tp_hybrid',
                    'binary': binary,
                    'kekule': kekule,
                    'peptide': {'_id': peptide['_id'],
                                'monomers': peptide['monomers']},
                    'template': template['_id']}
            self.result_data.append(doc)

    @staticmethod
    def create_tp_hybrid(peptide, template, logger=None):
        """
        Connects the peptide to the template through amide linkage with any primary amine on the peptide. This means a
        with a lysine will form two tp_hybrids, where one is connected through the N-terminal amine and the amine on the
        lysine.

        Args:
            template (rdkit Mol): The template molecule.
            peptide (rdkit Mol): The peptide molecule.
            n_term (str): The type of backbone the N-terminal monomer is made of.

        Returns:
            dict: A dictionary containing the tp_hybrid SMILES strings as keys and the kekule SMILES strings as values.
        """

        pep_mol = Chem.Mol(peptide['binary'])
        temp_mol = Chem.Mol(template['binary'])

        # find any primary amine or proline n-terminus
        matches = pep_mol.GetSubstructMatches(Chem.MolFromSmarts('[$([NH2]),$([NH;R]);!$([NH2]C(=O)*)]'))

        # for each eligible nitrogen form a connection
        results = {}
        for pairs in matches:

            # assign atom map number
            for atom_idx in pairs:
                atom = pep_mol.GetAtomWithIdx(atom_idx)
                if atom.GetSymbol() == 'N':  # primary amines only
                    atom.SetAtomMapNum(config.PEP_NITROGEN_MAP_NUM)
                    break

            # combine and record results
            try:
                tp_hybrid = Base.merge(temp_mol, pep_mol)
            except ValueError:
                if logger:
                    logger.exception(f'Sanitize Error! template = {Chem.MolToSmiles(temp_mol)}, '
                                      f'peptide = {Chem.MolToSmiles(pep_mol)}')
            else:
                atom.SetAtomMapNum(0)
                binary = tp_hybrid.ToBinary()
                Chem.Kekulize(tp_hybrid)
                results[binary] = Chem.MolToSmiles(tp_hybrid, kekuleSmiles=True)

        return NewTPHybrid(results, peptide, template)


########################################################################################################################
########################################################################################################################
########################################################################################################################

NewMacrocycles = namedtuple('NewMacrocycles', 'macrocycles tp_hybrid reaction')

class MacrocycleGenerator(Base):
    """
    Class for enumerating the candidate macrocycles by applying reactions to the tp_hybrids. Inherits from Base.

    Attributes:
        tp_hybrids (list): Contains the different side chains for generating regioisomers.
    """

    _defaults = config.DEFAULTS['MacrocycleGenerator']

    def __init__(self, logger=LOGGER, make_db_connection=True):
        """
        Initializer.

        Args:
        """

        # I/O
        super().__init__(logger, make_db_connection=True)

        # data
        self.tp_hybrids = []
        self.reactions = []

    def save_data(self):

        params = self._defaults['outputs']
        return self.to_mongo(params['col_macrocycles'])

    def load_data(self):
        """
        Overloaded method for loading input data.

        Returns:
            bool: True if successful.
        """

        try:
            params = self._defaults['inputs']
            self.tp_hybrids = self.from_mongo(params['col_tp_hyrbids'], {'type': 'tp_hybrid'})
            self.reactions = list(self.from_mongo(params['col_reactions'], {'type': {'$ne': 'template'}}))
        except Exception:
            self.logger.exception('Unexpected exception occured.')
        else:
            return True

        return False

    def generate(self):

        try:
            parent_sc_filter = partial(MacrocycleGenerator.is_applicable, attr='parent_side_chain')
            conn_atom_filter = partial(MacrocycleGenerator.is_applicable, attr='conn_atom_idx')
            template_filter = partial(MacrocycleGenerator.is_applicable, attr='template')
            first_filter = filter(parent_sc_filter, product(self.tp_hybrids, self.reactions))
            second_filter = filter(conn_atom_filter, first_filter)
            args = filter(template_filter, second_filter)

            with multiprocessing.Pool() as pool:
                results = pool.starmap_async(MacrocycleGenerator.apply_reaction, args)

                for macrocycles, tp_hybrid, reaction in results.get():
                    self.accumulate_data(macrocycles, tp_hybrid, reaction)
        except Exception:
            raise
        else:
            return True

        return False

    def generate_from_ids(self, tp_hybrid_ids, reaction_ids):
        try:
            params = self._defaults['inputs']
            self.tp_hybrids = self.from_mongo(params['col_tp_hyrbids'], {'type': 'tp_hybrid', '_id': {'$in': tp_hybrid_ids}})
            self.reactions = self.from_mongo(params['col_reactions'], {'type': 'reaction', '_id': {'$in': reaction_ids}})

            for tp_hybrid, reaction in product(self.tp_hybrids, self.reactions):
                if tp_hybrid['template'] == reaction['template']:
                    macrocycles, _, _ = MacrocycleGenerator.apply_reaction(tp_hybrid, reaction)
                    self.accumulate_data(macrocycles, tp_hybrid, reaction)
        except Exception:
            raise
        else:
            return True

        return False


    def accumulate_data(self, macrocycles, tp_hybrid, reaction):

        sc_id = reaction['side_chain']['_id']
        rxn_idx = reaction['rxn_atom_idx']
        for i, (binary, kekule) in enumerate(macrocycles.items()):
            doc = {
                '_id': tp_hybrid['_id'] + str(sc_id) + str(rxn_idx) + str(i),
                'type': 'macrocycle',
                'binary': binary,
                'kekule': kekule,
                'tp_hybrid': tp_hybrid['_id'],
                'reaction': {'_id': reaction['_id'],
                             'side_chain': sc_id,
                             'rxn_atom_idx': rxn_idx}
            }
            self.result_data.append(doc)

    @staticmethod
    def apply_reaction(tp_hybrid, reaction):
        """
        Applys all reaction templates to a reactant and collects all unique products

        Args:
            reactant (str): The reactant SMILES string
            rxn_docs (pymongo cursor): The collection of reaction documents containing reaction template SMARTS strings

        Returns:
            list: A list of SMARTS strings of all unique products
            list: A list of side chains that reacted in order to form the products
            list: A list of reacting atom indices in the side_chains that formed products
        """

        reactant = Chem.Mol(tp_hybrid['binary'])
        rxn = AllChem.ChemicalReaction(reaction['binary'])

        results = {}
        for products in rxn.RunReactants((reactant,)):
            for macrocycle in products:
                Chem.SanitizeMol(macrocycle)
                binary = macrocycle.ToBinary()
                Chem.Kekulize(macrocycle)
                results[binary] = Chem.MolToSmiles(macrocycle, kekuleSmiles=True)

        return NewMacrocycles(results, tp_hybrid, reaction)

    @staticmethod
    def is_applicable(args, attr):

        tp_hybrid, reaction = args
        if attr in ('parent_side_chain', 'conn_atom_idx'):
            for monomer in tp_hybrid['peptide']['monomers']:
                side_chain = monomer['side_chain']
                if side_chain is None or isinstance(side_chain, str):
                    continue
                if side_chain[attr] == reaction['side_chain'][attr]:
                    return True
        else:
            if tp_hybrid[attr] == reaction[attr]:
                return True

        return False

########################################################################################################################
########################################################################################################################
########################################################################################################################

NewConformers = namedtuple('NewConformers', 'macrocycle_doc macrocycle_mol energies rms ring_rms')
CleavedInfo = namedtuple('CleavedInfo', 'cleaved_atom1 cleaved_atom2 double_bond_flag')

class ConformerGenerator(Base):

    _defaults = config.DEFAULTS['ConformerGenerator']

    def __init__(self, logger=LOGGER, make_db_connection=True):

        # I/O
        super().__init__(logger, make_db_connection)

        # data
        self.macrocycles = []

    def save_data(self, mongo=True, pdb=False):

        try:
            params = self._defaults['outputs']
            for conformer in self.result_data:
                _id = conformer['_id']
                if pdb:
                    Chem.MolToPDBFile(Chem.Mol(conformer['binary']), os.path.join(params['fp_conformers'], _id + '.pdb'))
                if mongo:
                    del conformer['_id']
                    update = {'$set': {'conformer': conformer}}
                    self.mongo_db[params['col_conformers']].find_one_and_update({'_id': _id}, update)
        except Exception:
            raise
        else:
            return True

        return False

    def load_data(self):

        try:
            params = self._defaults['inputs']
            self.macrocycles = list(self.from_mongo(params['col_macrocycles'], {'type': 'macrocycle'})) #, 'regio_filter': True})
        except Exception:
            raise
        else:
            return True

        return False

    def generate(self, num_confs=1000, rms_threshold=1, seed=-1, max_iters=1000, d_min=0.35):
        try:
            embed_func = partial(ConformerGenerator.embed_conformers, num_confs=num_confs, rms_threshold=rms_threshold, seed=seed)
            optimize_func = partial(ConformerGenerator.optimize_conformers, max_iters=max_iters, d_min=d_min)
            with multiprocessing.Pool() as pool:
                embed_result = pool.map_async(embed_func, self.macrocycles)
                optimize_result = pool.starmap_async(optimize_func, embed_result.get())
                for macrocycle, conf_ids, energy, rms, _id in optimize_result.get():
                    self.accumulate_data(macrocycle, conf_ids, energy, rms, _id)
        except Exception:
            raise
        else:
            return True

        return False

    def generate_serial(self, num_confs=10, rms_threshold=0.5, seed=-1, max_iters=1000, d_min=0.35):

        try:
            for macrocycle in self.macrocycles:
                return ConformerGenerator.conformation_search(macrocycle, energy_diff=5)
                # tup = ConformerGenerator.embed_conformers(macrocycle, num_confs, seed, rms_threshold)
                # result = ConformerGenerator.optimize_conformers(*tup, max_iters=max_iters, d_min=d_min)
                # self.accumulate_data(result.macrocycle_doc, result.macrocycle_mol, result.conf_ids, result.energy, result.rms)
        except Exception:
            raise
        else:
            return True

        return False

    def accumulate_data(self, macrocycle, conf_ids, energy, rms, _id):

        doc = {
            '_id': _id,
            'binary': macrocycle.ToBinary(),
            'num_confs': num_confs,
            'energies': energies,
            'rms': rms,
            'ring_rms': ring_rms,
            'avg_energy': round(np.average(energies), 3),
            'avg_rms': round(np.average(rms), 3),
            'avg_ring_rms': round(np.average(ring_rms), 3)
        }
        self.result_data.append(doc)

    @staticmethod
    def conformation_search(macrocycle, repeats=5, num_confs_genetic=50, num_confs_keep=5, force_field='MMFF94s',
                            score='energy', min_rmsd=0.5, energy_diff=5, max_iters=1000, ring_size=10, granularity=5,
                            clash_threshold=0.9, distance_interval=[1.0, 2.5], seed=-1):

        storage_mol = Chem.Mol(macrocycle['binary'])


        old_smiles = Chem.MolToSmiles(storage_mol)
        cleavable_bonds = ConformerGenerator.get_cleavable_bonds(storage_mol, ring_size)
        dihedrals = ConformerGenerator.get_dihedral_atoms(storage_mol, ring_size)
        ring_atoms = ConformerGenerator.get_ring_atoms(storage_mol, ring_size)
        storage_mol = Chem.AddHs(storage_mol)

        # for each cleavable bond, perform algorithm
        opt_energies = {}
        min_energy = None
        params = AllChem.ETKDGv2()
        params.numThreads = 0
        for bond in cleavable_bonds:

            # cleave the bond and update the dihedral list
            linear_mol, cleaved_info = ConformerGenerator.cleave_bond(Chem.Mol(macrocycle['binary']), bond)
            linear_mol = Chem.AddHs(linear_mol)
            new_dihedrals = ConformerGenerator.update_dihedrals(linear_mol, dihedrals, cleaved_info.cleaved_atom1,
                                                                cleaved_info.cleaved_atom2)

            # use genetic algorithm to generate linear rotamers and optimize via force field then via dihedral rotations
            # and keep best results then repeat
            opt_linear_rotamers = []
            for i in range(repeats):
                rotamers = deepcopy(linear_mol)
                params.randomSeed = int(time()) if seed == -1 else seed
                while AllChem.EmbedMolecule(rotamers, params=params) < 0: params.randomSeed = int(time())
                ConformerGenerator.optimize_conformers(rotamers, force_field, max_iters)
                rotamers = ConformerGenerator.genetic_algorithm(rotamers, num_confs=num_confs_genetic, score=score)
                energies = ConformerGenerator.optimize_conformers(rotamers, force_field, max_iters)
                opt_linear_rotamers.extend(ConformerGenerator.optimize_linear_rotamers(rotamers,
                                           int(np.argmin(energies)), new_dihedrals, cleaved_info.cleaved_atom1,
                                           cleaved_info.cleaved_atom2, num_confs_keep, granularity, min_rmsd,
                                           clash_threshold, distance_interval, max_iters))

            # add best resulting rotamers to mol
            for optimized_linear in opt_linear_rotamers:
                linear_mol.AddConformer(optimized_linear, assignId=True)

            # reform bond and make sure the molecule hasnt changed
            macro_mol = ConformerGenerator.remake_bond(linear_mol, cleaved_info)
            if old_smiles != Chem.MolToSmiles(macro_mol):
                continue
            macro_mol = Chem.AddHs(macro_mol, addCoords=True)

            try:
                # optimize macrocycle and filter out conformers
                energies = ConformerGenerator.optimize_conformers(macro_mol, force_field, max_iters)
                ConformerGenerator.filter_conformers(macro_mol, energies, min_energy, energy_diff)
                mols = [ConformerGenerator.genetic_algorithm(macro_mol, conf_id=i, num_confs=num_confs_genetic,
                                                             score=score) for i in range(macro_mol.GetNumConformers())]
                macro_mol = ConformerGenerator.aggregate_conformers(mols)
                energies = ConformerGenerator.optimize_conformers(macro_mol, force_field, max_iters)
                ConformerGenerator.filter_conformers(macro_mol, energies, min_energy, energy_diff)

                # compare newly generated conformers to optimum conformers and if it is valid then add it to the list of
                # optimum conformers
                ConformerGenerator.evaluate_conformers(macro_mol, energies, storage_mol, opt_energies,
                                                       min_rmsd, max_iters)
                min_energy = min(opt_energies.values())
            except IndexError:  # number of conformers after filtering is 0
                continue

        # add conformers to opt_macrocycle in order of increasing energy
        energies, rms, ring_rms = [], [], []
        opt_macrocycle = Chem.AddHs(Chem.Mol(macrocycle['binary']))
        for conf_id, energy in sorted(opt_energies.items(), key=lambda x: x[1]):
            opt_macrocycle.AddConformer(storage_mol.GetConformer(conf_id), assignId=True)
            energies.append(energy)

        # align conformers
        AllChem.AlignMolConformers(opt_macrocycle, maxIters=max_iters, RMSlist=rms)
        AllChem.AlignMolConformers(opt_macrocycle, maxIters=max_iters, atomIds=ring_atoms, RMSlist=ring_rms)

        # remove temporary files
        ConformerGenerator.cleanup()

        return NewConformers(macrocycle, opt_macrocycle, energies, rms, ring_rms)

    @staticmethod
    def get_cleavable_bonds(macrocycle, ring_size=10):

        cleavable_bonds = []

        # identify the macrocycle rings' atoms and bonds
        macro_bonds = [ring for ring in macrocycle.GetRingInfo().BondRings() if len(ring) >= ring_size - 1] # bonds in ring = ring_size - 1
        macro_bonds = set().union(*macro_bonds)

        # identify chiral atoms
        chiral_atoms = [idx for idx, stereo in Chem.FindMolChiralCenters(macrocycle)]

        # find cleavable bonds
        for bond in macro_bonds:
            bond = macrocycle.GetBondWithIdx(bond)
            begin_atom = bond.GetBeginAtomIdx()
            end_atom = bond.GetEndAtomIdx()
            if bond.GetBondType() == Chem.BondType.SINGLE and begin_atom not in chiral_atoms \
               and end_atom not in chiral_atoms:
                cleavable_bonds.append(bond)

        return cleavable_bonds

    @staticmethod
    def get_dihedral_atoms(macrocycle, ring_size=10):

        dihedrals = {'cleaved_and_Hs': [],
                     'cleaved': [],
                     'other': []
                    }

        # get bonds in largest macrocycle ring
        macro_bonds = [list(ring) for ring in macrocycle.GetRingInfo().BondRings() if len(ring) >= ring_size - 1] # bonds in ring = ring_size - 1
        macro_bonds = sorted(macro_bonds, key=len, reverse=True)[0]
        macro_bonds.extend(macro_bonds[:2]) # duplicate first two bonds in list to the end; defines all possible dihedrals

        # find dihedrals - defined by bond patterns 1-2-3?4 or 1?2-3-4 where ? is any bond type
        for bond1, bond2, bond3 in window(macro_bonds, 3):
            bond1 = macrocycle.GetBondWithIdx(bond1)
            bond2 = macrocycle.GetBondWithIdx(bond2)
            bond3 = macrocycle.GetBondWithIdx(bond3)
            if bond2.GetBondType() == Chem.BondType.SINGLE and (bond1.GetBondType() == Chem.BondType.SINGLE or bond3.GetBondType() == Chem.BondType.SINGLE):

                # get correct ordering of dihedral atoms
                bond1_begin, bond1_end = bond1.GetBeginAtom(), bond1.GetEndAtom()
                bond3_begin, bond3_end = bond3.GetBeginAtom(), bond3.GetEndAtom()
                if bond1_begin.GetIdx() in [neighbor.GetIdx() for neighbor in bond3_begin.GetNeighbors()]:
                    dihedral = [bond1_end.GetIdx(), bond1_begin.GetIdx(), bond3_begin.GetIdx(), bond3_end.GetIdx()]
                elif bond1_begin.GetIdx() in [neighbor.GetIdx() for neighbor in bond3_end.GetNeighbors()]:
                    dihedral = [bond1_end.GetIdx(), bond1_begin.GetIdx(), bond3_end.GetIdx(), bond3_begin.GetIdx()]
                elif bond1_end.GetIdx() in [neighbor.GetIdx() for neighbor in bond3_begin.GetNeighbors()]:
                    dihedral = [bond1_begin.GetIdx(), bond1_end.GetIdx(), bond3_begin.GetIdx(), bond3_end.GetIdx()]
                else:
                    dihedral = [bond1_begin.GetIdx(), bond1_end.GetIdx(), bond3_end.GetIdx(), bond3_begin.GetIdx()]
                dihedrals['other'].append(dihedral)

        return dihedrals

    @staticmethod
    def get_ring_atoms(macrocycle, ring_size=10):
        ring_atoms = [ring for ring in macrocycle.GetRingInfo().AtomRings() if len(ring) >= ring_size - 1] # bonds in ring = ring_size - 1
        return  list(set().union(*ring_atoms))

    @staticmethod
    def cleave_bond(mol, bond):

        mol = Chem.RWMol(mol)

        # atom assignment must be done after conversion to RWMol
        cleaved_atom1 = mol.GetAtomWithIdx(bond.GetBeginAtomIdx())
        cleaved_atom2 = mol.GetAtomWithIdx(bond.GetEndAtomIdx())

        mol.RemoveBond(cleaved_atom1.GetIdx(), cleaved_atom2.GetIdx())

        # if alkene is adjacent to cleaved bond (cinnamoyl cation double bond)
        double_bond_flag = False
        for bond in cleaved_atom1.GetBonds() + cleaved_atom2.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE and not bond.IsInRing() and bond.GetBeginAtom().GetSymbol() == 'C' and bond.GetEndAtom().GetSymbol() == 'C':
                double_bond_flag = True
                break


        # adjust hydrogen counts
        cleaved_atom1.SetNumExplicitHs(1 + cleaved_atom1.GetTotalNumHs())
        cleaved_atom2.SetNumExplicitHs(1 + cleaved_atom2.GetTotalNumHs())
        Chem.SanitizeMol(mol)

        return mol, CleavedInfo(cleaved_atom1.GetIdx(), cleaved_atom2.GetIdx(), double_bond_flag)

    @staticmethod
    def remake_bond(mol, cleaved_info):

        # unpack values
        cleaved_atom1, cleaved_atom2 = cleaved_info.cleaved_atom1, cleaved_info.cleaved_atom2
        double_bond_flag = cleaved_info.double_bond_flag

        mol = Chem.RWMol(Chem.RemoveHs(mol, updateExplicitCount=True))

        cleaved_atom1 = mol.GetAtomWithIdx(cleaved_atom1)
        cleaved_atom2 = mol.GetAtomWithIdx(cleaved_atom2)

        # reset hydrogen counts
        cleaved_atom1.SetNumExplicitHs(cleaved_atom1.GetTotalNumHs() - 1)
        cleaved_atom2.SetNumExplicitHs(cleaved_atom2.GetTotalNumHs() - 1)

        mol.AddBond(cleaved_atom1.GetIdx(), cleaved_atom2.GetIdx(), Chem.BondType.SINGLE)

        # reset bond stereochemistry
        if double_bond_flag:
            new_bond = mol.GetBondBetweenAtoms(cleaved_atom1.GetIdx(), cleaved_atom2.GetIdx())
            new_bond.SetBondDir(Chem.BondDir.ENDUPRIGHT) # specific to cinnamoyl cation; causes double bond to be trans
            Chem.AssignStereochemistry(mol, force=True)

        Chem.SanitizeMol(mol)

        return mol

    @staticmethod
    def update_dihedrals(linear_mol, dihedrals, cleaved_atom1, cleaved_atom2):

        new_dihedrals = deepcopy(dihedrals)

        # find dihedrals that contain the bond that was cleaved and split that dihedral into two new lists that replaces
        # the opposite cleaved atom with one of the current cleaved atom's hydrogens
        for dihedral in dihedrals['other']:

            # both cleaved atoms are in dihedral
            if cleaved_atom1 in dihedral and cleaved_atom2 in dihedral:
                new_dh = deepcopy(dihedral)

                # cleaved_atom1 is a left terminal atom and cleaved_atom2 is a left central atom: ca1-ca2-x-x
                if cleaved_atom1 == dihedral[0]:
                    new_dh.remove(cleaved_atom1)
                    for neighbor in linear_mol.GetAtomWithIdx(cleaved_atom2).GetNeighbors():
                        if neighbor.GetSymbol() == 'H':
                            new_dh.insert(0, neighbor.GetIdx())
                            break
                    new_dihedrals['cleaved_and_Hs'].append(new_dh)

                # cleaved_atom1 is a right terminal atom and cleaved_atom2 is a right central atom: x-x-ca2-ca1
                elif cleaved_atom1 == dihedral[-1]:
                    new_dh.remove(cleaved_atom1)
                    for neighbor in linear_mol.GetAtomWithIdx(cleaved_atom2).GetNeighbors():
                        if neighbor.GetSymbol() == 'H':
                            new_dh.append(neighbor.GetIdx())
                            break
                    new_dihedrals['cleaved_and_Hs'].append(new_dh)

                # cleaved_atom2 is a left terminal atom and cleaved_atom1 is a left central atom: ca2-ca1-x-x
                elif cleaved_atom2 == dihedral[0]:
                    new_dh.remove(cleaved_atom2)
                    for neighbor in linear_mol.GetAtomWithIdx(cleaved_atom1).GetNeighbors():
                        if neighbor.GetSymbol() == 'H':
                            new_dh.insert(0, neighbor.GetIdx())
                            break
                    new_dihedrals['cleaved_and_Hs'].append(new_dh)

                # cleaved_atom2 is a right terminal atom and cleaved_atom1 is a right central atom: x-x-ca1-ca2
                elif cleaved_atom2 == dihedral[-1]:
                    new_dh.remove(cleaved_atom2)
                    for neighbor in linear_mol.GetAtomWithIdx(cleaved_atom1).GetNeighbors():
                        if neighbor.GetSymbol() == 'H':
                            new_dh.append(neighbor.GetIdx())
                            break
                    new_dihedrals['cleaved_and_Hs'].append(new_dh)

                # nothing special to be done for when cleaved_atom1 and cleaved_atom2 are both central: x-ca1-ca2-x or
                # x-ca2-ca1-x

                # remove dihedral from original list
                new_dihedrals['other'].remove(dihedral)

            # only one cleaved atom in dihedral
            elif cleaved_atom1 in dihedral or cleaved_atom2 in dihedral:
                new_dihedrals['cleaved'].append(dihedral)
                new_dihedrals['other'].remove(dihedral)


        return new_dihedrals

    @staticmethod
    def optimize_conformers(mol, force_field='MMFF94s', max_iters=1000):

        mol_props = AllChem.MMFFGetMoleculeProperties(mol)
        mol_props.SetMMFFVariant(force_field)
        force_fields = list(map(lambda x: AllChem.MMFFGetMoleculeForceField(mol, mol_props, confId=x, ignoreInterfragInteractions=False), range(mol.GetNumConformers())))

        convergence, energy = zip(*AllChem.MMFFOptimizeMoleculeConfs(mol, mmffVariant=force_field, maxIters=max_iters, numThreads=0, ignoreInterfragInteractions=False))
        convergence = list(convergence)
        energy = list(energy)
        for non_converged_id in np.flatnonzero(convergence):
            while not convergence[non_converged_id]:
                convergence[non_converged_id] = force_field[non_converged_id].Minimize(50)
                energy[non_converged_id] = force_field[non_converged_id].CalcEnergy()

        return list(map(lambda x: x.CalcEnergy(), force_fields))

    @staticmethod
    def genetic_algorithm(mol, conf_id=0, num_confs=50, score='energy', remove_Hs=False):

        outputs = ConformerGenerator._defaults['outputs']
        mol_file, results_file = outputs['tmp_molecule'], outputs['tmp_genetic_results']

        command = f'obabel {mol_file} -O {results_file} --conformer --nconf {num_confs} --score {score} \
                    --writeconformers &> /dev/null'
        write_mol(mol, mol_file, conf_id=conf_id)
        os.system(command)
        return ConformerGenerator.aggregate_conformers([mol for mol in Chem.SDMolSupplier(results_file,
                                                                                          removeHs=remove_Hs)])

    @staticmethod
    def aggregate_conformers(mols):

        # get mol to add rest of conformers to
        mol = mols.pop()

        # add conformers from rest of mols to mol
        for conformer in [conf for m in mols for conf in m.GetConformers()]:
            mol.AddConformer(conformer, assignId=True)

        return mol

    @staticmethod
    def optimize_linear_rotamers(linear_mol, conf_id, dihedrals, cleaved_atom1, cleaved_atom2, num_confs=5,
                                 granularity=5, min_rmsd=0.5, clash_threshold=0.9, distance_interval=[1.0, 2.5],
                                 max_iters=1000):

        mast_mol = deepcopy(linear_mol)
        mast_mol.RemoveAllConformers()
        optimized_linear_confs, distances = [], []
        linear_conf = linear_mol.GetConformer(conf_id)

        # generate length 2 combinations for dihedrals that don't contain cleaved atoms and get the resulting
        # distances between the two cleaved atoms after applying various angles to those dihedrals. Sort the results
        # based on distance
        for dihedral1, dihedral2 in combinations(dihedrals['other'], 2):
            ini_dihedral1 = AllChem.GetDihedralDeg(linear_conf, dihedral1[0], dihedral1[1], dihedral1[2], dihedral1[3])
            ini_dihedral2 = AllChem.GetDihedralDeg(linear_conf, dihedral2[0], dihedral2[1], dihedral2[2], dihedral2[3])
            dist = ConformerGenerator.get_distance(linear_conf, cleaved_atom1, cleaved_atom2)
            if distance_interval[0] < dist < distance_interval[1]:
                distances.append([dist, ini_dihedral1, dihedral1, ini_dihedral2, dihedral2])

            angle1, angle2 = 0, 0
            while angle1 < 360:
                AllChem.SetDihedralDeg(linear_conf, dihedral1[0], dihedral1[1], dihedral1[2], dihedral1[3], angle1)
                while angle2 < 360:
                    AllChem.SetDihedralDeg(linear_conf, dihedral2[0], dihedral2[1], dihedral2[2], dihedral2[3], angle2)
                    dist = ConformerGenerator.get_distance(linear_conf, cleaved_atom1, cleaved_atom2)
                    if distance_interval[0] < dist < distance_interval[1]:
                        distances.append([dist, angle1, dihedral1, angle2, dihedral2])
                    angle2 += granularity
                angle1 += granularity

            # reset dihedrals
            AllChem.SetDihedralDeg(linear_conf, dihedral1[0], dihedral1[1], dihedral1[2], dihedral1[3], ini_dihedral1)
            AllChem.SetDihedralDeg(linear_conf, dihedral2[0], dihedral2[1], dihedral2[2], dihedral2[3], ini_dihedral2)
        distances.sort(key=lambda x: x[0])

        # starting with the dihedral combinations that minimized the distance between cleaved atoms the most, find
        # the optimimum angles for dihedrals that contain cleaved atoms and no hydrogens, then for dihedrals that
        # contain cleaved atoms and hydrogens, until desired number of conformers has been generated
        for distance in distances:
            linear_mol_copy = deepcopy(linear_mol)
            linear_conf = linear_mol_copy.GetConformer(conf_id)

            # set starting dihedrals
            AllChem.SetDihedralDeg(linear_conf, distance[2][0], distance[2][1], distance[2][2], distance[2][3], distance[1])
            AllChem.SetDihedralDeg(linear_conf, distance[4][0], distance[4][1], distance[4][2], distance[4][3], distance[3])

            # if no clashes are detected optimize continue optimization
            matrix = Chem.Get3DDistanceMatrix(linear_mol, confId=conf_id).flatten()
            matrix = matrix[matrix > 0]
            if sum(matrix < clash_threshold) == 0:

                # optimize dihedrals
                ConformerGenerator.optimize_dihedrals(linear_conf, dihedrals['cleaved'], cleaved_atom1, cleaved_atom2, granularity)
                ConformerGenerator.optimize_dihedrals(linear_conf, dihedrals['cleaved_and_Hs'], cleaved_atom1, cleaved_atom2, granularity)

                for ref_conf in range(mast_mol.GetNumConformers()):
                    rms = AllChem.AlignMol(linear_mol_copy, mast_mol, conf_id, ref_conf, maxIters=max_iters)
                    if rms < min_rmsd:
                        break
                else:
                    optimized_linear_confs.append(linear_conf)
                    mast_mol.AddConformer(linear_conf, assignId=True)

                # return when num_confs valid conformers has been obtained
                if len(optimized_linear_confs) == num_confs:
                    break

        return optimized_linear_confs

    @staticmethod
    def optimize_dihedrals(conformer, dihedrals, cleaved_atom1, cleaved_atom2, granularity=5):

        for dihedral in dihedrals:
            best_dist = ConformerGenerator.get_distance(conformer, cleaved_atom1, cleaved_atom2)
            best_angle = AllChem.GetDihedralDeg(conformer, dihedral[0], dihedral[1], dihedral[2], dihedral[3])
            angle = 0
            while angle < 360:
                AllChem.SetDihedralDeg(conformer, dihedral[0], dihedral[1], dihedral[2], dihedral[3], angle)
                dist = ConformerGenerator.get_distance(conformer, cleaved_atom1, cleaved_atom2)
                if dist < best_dist:
                    best_dist = dist
                    best_angle = angle
                angle += granularity
            AllChem.SetDihedralDeg(conformer, dihedral[0], dihedral[1], dihedral[2], dihedral[3], best_angle)

    @staticmethod
    def get_distance(mol, atom1, atom2):

        atom1_position = mol.GetAtomPosition(atom1)
        atom2_position = mol.GetAtomPosition(atom2)
        return atom1_position.Distance(atom2_position)

    @staticmethod
    def evaluate_conformers(mol, energies, optimal_mol, opt_energies, min_rmsd=0.5, max_iters=1000):

        for i, macro_conf in enumerate(mol.GetConformers()):
            similar_confs = []
            for opt_conf in optimal_mol.GetConformers():
                rmsd = AllChem.AlignMol(mol, optimal_mol, macro_conf.GetId(), opt_conf.GetId(), maxIters=max_iters)
                if rmsd < min_rmsd:
                    similar_confs.append(opt_conf.GetId())
            similar_energies = [opt_energies[conf_id] for conf_id in similar_confs]
            similar_energies.append(energies[i])
            try:
                max_id = similar_confs[np.argmax(similar_energies)]
                optimal_mol.RemoveConformer(max_id)
                conf_id = optimal_mol.AddConformer(macro_conf, assignId=True)
                optimal_mol.GetConformer(conf_id).SetId(max_id)
                opt_energies[max_id] = energies[i]
            except IndexError:
                if len(similar_confs) == 0:
                    conf_id = optimal_mol.AddConformer(macro_conf, assignId=True)
                    opt_energies[conf_id] = energies[i]

    @staticmethod
    def filter_conformers(mol, energies, min_energy=None, energy_diff=5):

        remove_flag = False
        if min_energy is None:
            min_energy = min(energies)

        # filter high energy confs
        for conf_id, energy in zip(range(mol.GetNumConformers()), deepcopy(energies)):
            if energy > min_energy + energy_diff:
                mol.RemoveConformer(conf_id)
                energies.remove(energy)
                remove_flag = True

        # reset conf_ids if conformers have been filtered out
        if remove_flag:
            for i, conformer in enumerate(mol.GetConformers()):
                conformer.SetId(i)

    @staticmethod
    def cleanup():

        outputs = ConformerGenerator._defaults['outputs']
        for file in (outputs['tmp_molecule'], outputs['tmp_genetic_results']):
            if os.path.exists(file):
                os.remove(file)

