import macrocycles.utils as utils
from logging import INFO
import macrocycles.config as config
import multiprocessing
from itertools import product
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem
import csv
from collections import defaultdict
import json
from bson import json_util

LOGGER = utils.create_logger(__name__, INFO)


class RegioSQMFilter(utils.Base):

    _defaults = config.DEFAULTS['RegioSQMFilter']

    def __init__(self, logger=LOGGER, make_db_connection=True):

        super().__init__(logger, make_db_connection)

        self.macrocycles = []
        self.predictions = []

    def create_smiles_file(self, filepath):

        mols = self.mongo_db[config.COL1].find({'type': 'side_chain', 'connection': 'methyl'})
        with open(filepath, 'w') as file:
            for i, mol in enumerate(mols):
                file.write('sc' + str(i) + '\t' + mol['kekule'] + '\n')

    def import_predictions(self, smiles_filepath, csv_filepath, solvent, cut_off):
        smiles_dict = {}
        with open(smiles_filepath, 'r') as file:
            for line in file.readlines():
                line = line.split('\t')
                smiles_dict[line[0]] = line[1].rstrip()  # line[0] = identifier, line[1] = SMILES string

        collection = []
        with open(csv_filepath, 'r') as file:
            data = csv.reader(file, delimiter=',')
            for i, row in enumerate(data):
                doc = {}
                if i % 3 == 0:
                    mol = Chem.MolFromSmiles(smiles_dict[row[0]])
                    Chem.Kekulize(mol)
                    kekule = Chem.MolToSmiles(mol, kekuleSmiles=True)
                elif i % 3 == 1:
                    prediction = [int(row[j]) for j in range(0, len(row), 2)]
                else:
                    side_chain = self.mongo_db[config.COL1].find_one({'kekule': kekule})
                    doc = {'_id': side_chain['_id'] + 'r',
                           'type': 'regiosqm',
                           'kekule': kekule,
                           'predictions': prediction,
                           'solvent': solvent,
                           'cut_off': cut_off,
                           'parent_side_chain': side_chain['parent_side_chain']['_id'],
                           'conn_atom_idx': side_chain['conn_atom_idx']}
                    collection.append(doc)

        self.mongo_db[config.COL3].insert_many(collection)

    def save_data(self, to_mongo=True, to_json=False):

        try:
            params = self._defaults['outputs']

            if to_mongo:
                for macrocycle in self.result_data:
                    update = {'$set': {'tractable': macrocycle['tractable']}}
                    self.mongo_db[params['col_filtered']].find_one_and_update({'_id': macrocycle['_id']}, update)

            if to_json:
                with open(params['json_filtered'], 'w') as file:
                    json.dump(json.loads(json_util.dumps(self.macrocycles)), file)

        except Exception:
            raise
        else:
            return True

        return False

    def load_data(self, source='mongo', solvent='nitromethane', cut_off=3):

        try:
            params = self._defaults['inputs']
            if source == 'mongo':
                self.macrocycles = list(self.from_mongo(params['col_macrocycles'], {'type': 'macrocycle'}))
                self.predictions = list(self.from_mongo(
                    params['col_filters'], {'type': 'regiosqm', 'solvent': solvent, 'cut_off': cut_off}))
            elif source == 'json':
                self.macrocycles = self.from_json(params['json_macrocycles'])
                self.predictions = self.from_json(params['json_predictions'])
            else:
                self.logger.warning(f'The source type \'{source}\' for input data is unrecognized')
        except Exception:
            raise
        else:
            return True

        return False

    def filter_macrocycles(self):

        try:
            args = self.get_args()
            with multiprocessing.Pool() as pool:
                results = pool.starmap_async(RegioSQMFilter.apply_filter, args)

                for macrocycle in results.get():
                    self.result_data.append(macrocycle)

        except Exception:
            raise
        else:
            return True

        return False

    @staticmethod
    def apply_filter(macrocycle, filter_docs):

        flag = False
        for reaction in macrocycle['reactions']:
            for filter_doc in filter_docs[reaction['parent_side_chain']]:
                if reaction['conn_atom_idx'] == filter_doc['conn_atom_idx']:
                    macrocycle.update({'tractable': reaction['rxn_atom_idx'] in filter_doc['predictions']})
                    if not macrocycle['tractable']:
                        flag = True
                    break

            if flag:
                break

        return macrocycle

    def get_args(self):

        # hash predictions based on parent side chain _id
        predictions = defaultdict(list)
        for prediction in self.predictions:
            predictions[prediction['parent_side_chain']].append(prediction)

        # find applicable prediction documents for each macrocycle
        applicable_reactions = {'friedel_crafts', 'pictet_spangler', 'pyrrolo_indolene'}
        for macrocycle in self.macrocycles:
            reaction_types = set(reaction['type'] for reaction in macrocycle['reactions'])
            if reaction_types.intersection(applicable_reactions):
                applicable_predictions = {}
                parent_side_chain_ids = [reaction['parent_side_chain'] for reaction in macrocycle['reactions']]
                for parent_side_chain_id in parent_side_chain_ids:
                    applicable_predictions[parent_side_chain_id] = predictions[parent_side_chain_id]

                yield (macrocycle, applicable_predictions)


class pKaFilter(utils.Base):

    def __init__(self):

        self.heterocycles = []

    def import_predictions(self, filepath):

        data = []
        with open(filepath, 'r') as file:
            solvent = file.readline().strip('\n')
            for line in file.readlines():
                line = line.strip('\n').replace(' ', '')
                smiles, *vals = line.split(';')
                mol = Chem.MolFromSmiles(smiles)
                predictions = {}
                for atom in mol.GetAtoms():
                    atom_map_num = atom.GetAtomMapNum()
                    if atom_map_num:
                        avg, std = vals[atom_map_num - 1].split(',')
                        predictions[atom.GetIndex()] = {'avg': float(avg),
                                                        'std': float(std)}
                        atom.SetAtomMapNum(0)
                Chem.Kekulize(mol)
                data.append({'_id': None,
                             'type': 'pka',
                             'kekule': Chem.MolToSmiles(mol, kekuleSmiles=True),
                             'prediction': predictions,
                             'solvent': solvent})

        self.mongo_db['filters'].insert_many(data)

    def save_data(self):
        pass

    def load_data(self):
        pass


# def indole_hetero_atom_filter(collection):
#     """
#     For each reactant in candidates, apply a filter that selects the enumerated candidates that result from an indole
#     nitrogen reacting to close the macrocycle. This filter is generalized to any indole derivative. The filtered
#     candidates are stored in a dictionary with the following hierarchy:
#         - reacting side chain
#             - candidates

#     Args:
#         collection (iterable): An iterable object containing the candidate data, such as a pymongo cursor on the
#             candidates collection

#     Returns:
#         list: A list containing tuples of the reactant and the sorted dictionary of filtered candidates
#     """

#     results = []
#     for doc in collection:

#         # initialize dict based on unique side chains that reacted to close macrocycle
#         filtered_cands = {}
#         for reacting_sc in list(set(doc['reacting_side_chains'])):
#             filtered_cands[reacting_sc] = []

#         # apply filter
#         for product, side_chain in zip(doc['products'], doc['reacting_side_chains']):
#             product = Chem.MolFromSmiles(product)
#             if product.HasSubstructMatch(Chem.MolFromSmarts(INDOLE_SMARTS)):
#                 filtered_cands[side_chain].append(Chem.MolToSmiles(product))

#         # delete entries with no values
#         filtered_cands = clean_results(filtered_cands)

#         results.append((doc['reactant'], filtered_cands))

#     return results
