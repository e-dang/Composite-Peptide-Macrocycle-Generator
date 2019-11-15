import macrocycles.utils as utils
from logging import INFO
import macrocycles.config as config
import multiprocessing
from itertools import product
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem
import csv

LOGGER = utils.create_logger(__name__, INFO)


class RegioSQMFilter(utils.Base):

    _defaults = config.DEFAULTS['RegioSQMFilter']

    def __init__(self, logger=LOGGER, make_db_connection=True):

        super().__init__(logger, make_db_connection)

        self.macrocycles = []
        self.regio_predictions = []

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
                elif i % 3 == 2:
                    if len(row) != 0:
                        print('error')
                    doc = {'_id': self.mongo_db['molecules'].find_one({'kekule': kekule}, projection={'_id': 1})['_id'],
                           'type': 'regiosqm',
                           'kekule': kekule,
                           'prediction': prediction,
                           'solvent': solvent,
                           'cut_off': cut_off}
                    collection.append(doc)

        self.mongo_db['filters'].insert_many(collection)

    def save_data(self):

        try:
            params = self._defaults['outputs']
            for macrocycle in self.result_data:
                update = {'$set': {'regio_filter': macrocycle['regio_filter']}}
                self.mongo_db[params['col_filtered']].find_one_and_update({'_id': macrocycle['_id']}, update)
        except Exception:
            raise
        else:
            return True

        return False

        params = self._defaults['outputs']
        return self.to_mongo(params['col_regiosqm'])

    def load_data(self, solvent, cut_off):

        try:
            params = self._defaults['inputs']
            self.macrocycles = self.from_mongo(params['col_macrocycles'], {'type': 'macrocycle'})
            self.regio_predictions = self.from_mongo(
                params['col_filters'], {'type': 'regiosqm', 'solvent': solvent, 'cut_off': cut_off})
        except Exception:
            raise
        else:
            return True

        return False

    def filter_macrocycles(self):

        try:
            args = filter(RegioSQMFilter.is_applicable, product(self.macrocycles, self.regio_predictions))

            with multiprocessing.Pool() as pool:
                results = pool.starmap_async(RegioSQMFilter.apply_filters, args)

                for macrocycle, filter_doc in results.get():
                    self.result_data.append(macrocycle)

        except Exception:
            raise
        else:
            return True

        return False

    @staticmethod
    def apply_filters(macrocycle, filter_doc):

        try:
            if not macrocycle['regio_filter']:
                return macrocycle, filter_doc
        except KeyError:
            for reaction in macrocycle['reaction']:
                if reaction['side_chain'] == filter_doc['_id']:
                    macrocycle.update({'regio_filter': reaction['rxn_atom_idx'] in filter_doc['prediction']})
                    if not macrocycle['regio_filter']:
                        break

        return macrocycle, filter_doc

    @staticmethod
    def is_applicable(args):

        macrocycle, filter_doc = args
        for reaction in macrocycle['reaction']:
            if reaction['side_chain'] == filter_doc['_id']:
                return True

        return False


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
                        avg, std = vals[atom_map_num].split(',')
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
