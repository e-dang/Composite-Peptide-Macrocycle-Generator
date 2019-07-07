"""
Written by Eric Dang.
github: https://github.com/e-dang
email: edang830@gmail.com
"""

from logging import INFO
from collections import ChainMap, namedtuple
from copy import deepcopy
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import json
from bson import json_util
from itertools import cycle, chain, islice, product, dropwhile
from functools import partial
import multiprocessing
from random import sample
from pprint import pprint

from macrocycles.exceptions import MissingMapNumberError, AtomSearchError
from macrocycles.utils import Base, create_logger, read_mols, get_user_approval, get_user_atom_idx, atom_to_wildcard, ranges
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
        self.load_rtemplates()
        self.result_data = []
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

        if group is None:
            group = input(f'Enter group name for monomers in file {fp}: ')

        if collection is None:
            collection = self._defaults['outputs']['col_monomers']

        template_doc = {'_id': None,
                        'type': 'monomer',
                        'binary': None,
                        'kekule': None,
                        'backbone': backbone,
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
            collection = self._defaults['outputs']['col_templates']

        template_doc = {'_id': None,
                        'type': 'template',
                        'binary': None,
                        'kekule': None,
                        'original': None}

        if self.from_sdf(fp, template_doc, DataInitializer.set_template_properties):
            return self.to_mongo(collection, create_id=False)

        return False

    @staticmethod
    def set_template_properties(template, doc):

        # get SMILES of original template before modification
        doc['original'] = Chem.MolToSmiles(template)

        # get _id from user
        Draw.ShowMol(template)
        doc['_id'] = input('Enter the name of the template: ')

        # set atom map number for connecting atom between template and peptide
        question = 'Enter the index of the carbon atom that will form amide linkage to the peptides: '
        substruct = DataInitializer.generalize_substruct(template, config.TEMP_CARBON_MAP_NUM, question)
        template = Chem.DeleteSubstructs(template, Chem.MolFromSmiles(substruct))

        return template, doc

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
                side_chain = Base.merge(parent_mol, conn_mol, config.PSC_MAP_NUM, config.CONN_MAP_NUM)
            except ValueError:
                if logger:
                    logger.exception(f'Sanitize error! Parent Side Chain: {Chem.MolToSmiles(parent_mol)}')
                else:
                    print('error')
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
                monomer = Base.merge(sc_mol, bb_mol, config.SC_MAP_NUM, config.BB_MAP_NUM, stereo)
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
            monomer_prod = product(self.monomers, repeat=length)
            if num_peptides:
                if random:
                    monomers = sample(monomer_prod, num_peptides)
                else:
                    monomers = islice(monomer_prod, num_peptides)
            else:
                start, stop = ranges(len(self.monomers) ** length, num_jobs)[job_id - 1]
                monomers = islice(monomer_prod, start, stop)

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
                monomers = self.from_mongo(params['col_monomers'], {'type': 'monomer', '_id': {'$in': monomers}})
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
        # for monomer in monomers:
        #     monomer_data.append({'_id': monomer['_id'], 'side_chain': monomer['side_chain']})

        # monomer_ids = [(monomer['_id'], monomer['side_chain']['conn_atom_idx']) for monomer in monomers]
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
                peptide = Base.merge(peptide, monomer_mol, config.PEP_CARBON_MAP_NUM, config.MONO_NITROGEN_MAP_NUM)
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
                tp_hybrid = Base.merge(temp_mol, pep_mol, config.TEMP_CARBON_MAP_NUM, config.PEP_NITROGEN_MAP_NUM)
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
        return self.to_mongo(params['col_macrocycles'], {'type': 'macrocycle'})

    def load_data(self):
        """
        Overloaded method for loading input data.

        Returns:
            bool: True if successful.
        """

        try:
            params = self._defaults['inputs']
            self.tp_hybrids = self.from_mongo(params['col_tp_hyrbids'], {'type': 'tp_hybrid'})
            self.reactions = self.from_mongo(params['col_reactions'], {'type': {'$ne': 'template'}})
        except Exception:
            self.logger.exception('Unexpected exception occured.')
        else:
            return True

        return False

    # def generate_candidates(self):

    #     for tp_hybrid in self.tp_hybrids:
    #         reactant = tp_hybrid['smiles']
    #         print('', reactant)
    #         for i, monomer in enumerate(set(tp_hybrid['peptide']['monomers'])):
    #             print('\t', monomer)
    #             monomer = self.mongo_db['molecules'].find_one({'ID': monomer})
    #             # print(monomer)
    #             reactions = list(self.mongo_db['reactions'].find(
    #                 {'type': 'reaction', 'side_chain.side_chain.parent.ID': monomer['side_chain']['parent']['ID']}))
    #             reactions.extend(list(self.mongo_db['reactions'].find(
    #                 {'type': 'reaction', 'side_chain': monomer['side_chain']['ID']}
    #             )))
    #             # print(reactions[0])
    #             print(len(reactions))
    #             # print(reactions)
    #             products, reacting_side_chains, atom_idxs = self.apply_reaction(reactant, reactions)
    #             self.accumulate_data(products, reacting_side_chains, atom_idxs, tp_hybrid, reactions, i)

    #     return True

    def generate(self):
        for tp_hybrid in self.tp_hybrids:
            # lookups = defaultdict(list)
            # unique = ((monomer['side_chain']['parent_side_chain'], monomer['side_chain']['conn_atom_idx']) for monomer in tp_hybrid['peptide']['monomers'])
            # for key, val in unique:
            #     lookups[key].append(val)
            for monomer in tp_hybrid['peptide']['monomers']:
                side_chain = monomer['side_chain']
                if side_chain is not None:
                    reactions = self.from_mongo('reactions', {'side_chain.parent_side_chain': side_chain['parent_side_chain'], 'side_chain.conn_atom_idx': side_chain['conn_atom_idx']})
                    for reaction in reactions:
                        macrocycles, _, _ = MacrocycleGenerator.apply_reaction(tp_hybrid, reaction)
                        self.accumulate_data(macrocycles, tp_hybrid, reaction)
        return True

    def generate_parallel(self):

        try:
            args = self.pair_parameters()

            with multiprocessing.Pool() as pool:
                results = pool.starmap_async(MacrocycleGenerator.apply_reaction, args)

                for macrocycles, tp_hybrid, reaction in results.get():
                    self.accumulate_data(macrocycles, tp_hybrid, reaction)
        except Exception:
            raise
        else:
            return True

        return False


    def accumulate_data(self, macrocycles, tp_hybrid, reaction):

        for i, (binary, kekule) in enumerate(macrocycles.items()):
            doc = {
                '_id': tp_hybrid['_id'] + str(i),
                'type': 'macrocycle',
                'binary': binary,
                'kekule': kekule,
                'tp_hybrid': tp_hybrid['_id'],
                'reaction': {'_id': reaction['_id'],
                             'side_chain': reaction['side_chain']['_id'],
                             'rxn_atom_idx': reaction['rxn_atom_idx']}
            }
            self.result_data.append(doc)

    # def accumulate_data(self, macrocycles, reacting_side_chains, atom_idxs, tp_hybrid, reactions, i):

    #     # for j, (macrocycle, reacting_side_chain, atom_idx) in enumerate(zip(macrocycles, reacting_side_chains, atom_idxs)):
    #     doc = {
    #         'ID': 'c' + str(i) + tp_hybrid['ID'],
    #         'type': 'candidate',
    #         'smiles': macrocycles,
    #         'kekule': '',
    #         'reacting_side_chains': reacting_side_chains,
    #         'atom_idxs': atom_idxs,
    #         'tp_hybrid': tp_hybrid
    #         # 'reaction': reactions
    #     }
    #     self.result_data.append(doc)

    def pair_parameters(self):

        parent_side_chains = list(self.from_mongo('molecules', {'type': 'parent_side_chain'}))
        for parent_sc in parent_side_chains:
            tp_hybrids = dropwhile(lambda x: MacrocycleGenerator.has_parent_side_chain(x, parent_sc), self.tp_hybrids)
            reactions = dropwhile(lambda x: MacrocycleGenerator.has_parent_side_chain(x, parent_sc), self.reactions)
            for tp_hybrid in tp_hybrids:
                applicable_rxns = dropwhile(lambda x: MacrocycleGenerator.has_same_conn_idx(x, tp_hybrid), reactions)
                yield cycle(tp_hybrid), applicable_rxns


    @staticmethod
    def has_parent_side_chain(tp_hybrid, parent_sc):

        for monomer in tp_hybrid['peptide']['monomers']:
            if monomer['side_chain']['parent_side_chain'] == parent_sc['_id']:
                return True

        return False

    @staticmethod
    def has_same_conn_idx(reactions, tp_hybrid):

        for reaction in reactions:
            if reaction['side_chain']['conn_atom_idx'] == tp_hybrid['side_chain']['conn_atom_idx']:
                return True

        return False

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

        products = rxn.RunReactants((reactant,))
        for product in products:
            Chem.SanitizeMol(product[0])
            binary = product.ToBinary()
            Chem.Kekulize(product)
            results[binary] = Chem.MolToSmiles(product, kekuleSmiles=True)

        return NewMacrocycles(results, tp_hybrid, reaction)

        # # apply reactions
        # reaction_info = {}
        # for doc in reactions:
        #     rxn = AllChem.ReactionFromSmarts(doc['binary'])
        #     prod = rxn.RunReactants((reactant,))
        #     # print(prod)

        #     # get unique products
        #     for mol in prod:
        #         Chem.SanitizeMol(mol[0])
        #         smiles = Chem.MolToSmiles(mol[0])
        #         # print(smiles)
        #         if smiles not in reaction_info.keys():
        #             try:
        #                 reaction_info[smiles] = (doc['side_chain'], doc['side_chain']['rxn_atom_idx'])
        #             except:
        #                 reaction_info[smiles] = (doc['side_chain'], None)

        # # print(reaction_info)
        # unique_prods, sc_and_idx = zip(*[(key, value) for key, value in reaction_info.items()])
        # reacting_side_chains, atom_idxs = zip(*list(sc_and_idx))
        # # print(reacting_side_chains)
        # # print(atom_idxs)

        # return unique_prods, reacting_side_chains, atom_idxs