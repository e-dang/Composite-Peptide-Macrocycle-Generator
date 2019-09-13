import macrocycles.utils as utils
from logging import INFO
import macrocycles.config as config
import multiprocessing
from itertools import product

LOGGER = utils.create_logger(__name__, INFO)


class RegioSQMFilter(utils.Base):

    _defaults = config.DEFAULTS['RegioSQMFilter']

    def __init__(self, logger=LOGGER, make_db_connection=True):

        super().__init__(logger, make_db_connection)

        self.macrocycles = []
        self.regio_predictions = []

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
        # print(macrocycle['reaction'])
        # print(filter_doc)
        for reaction in macrocycle['reaction']:
            if reaction['side_chain'] == filter_doc['_id']:
                return True

        return False


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
