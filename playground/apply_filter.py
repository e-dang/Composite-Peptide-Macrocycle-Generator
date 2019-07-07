import argparse
from copy import deepcopy

from rdkit import Chem
from tqdm import tqdm
from pprint import pprint

from macrocycles.utils import MongoDataBase

REGIOSQM_FILTER = 'regioSQM'
INDOLE_FILTER = 'indole'

INDOLE_SMARTS = '*/C=C/Cn1aaa2aaaaa21'


# def regiosqm_filter(collection):
#     """
#     For each doc in the collection, the regioSQM filter is applied to the enumerated candidates resulting from each
#     unique side chain reacting to close the macrocycle. The filtered candidates are stored in a dictionary with the
#     following hierarchy:
#         - reacting side chain
#             - threshold
#                 - candidates
#                 - heats

#     Args:
#         collection (iterable): An iterable object containing the candidate data, such as a pymongo cursor on the
#             candidates collection

#     Returns:
#         list: A list containing tuples of the reactant and the sorted dictionary of filtered candidates
#     """

#     results = []
#     for doc in collection:

#         # get all unique side_chains that reacted in the reactant and corresponding regioSQM predictions
#         unique_sc = list(set(doc['reacting_side_chains']))
#         regio_filter = Database(db='rxn_templates').find('regioSQM', {'side_chain': {'$in': unique_sc}}, {'_id': 0})

#         # initialize thresholds for each unique side_chain
#         filtered_cands = {}
#         for reacting_sc in unique_sc:
#             filtered_cands[reacting_sc] = {}
#             filtered_cands[reacting_sc]['threshold_1'] = []
#             filtered_cands[reacting_sc]['threshold_2'] = []

#         # place enumerated products in appropriate theshold bin
#         for product, side_chain, atom_idx in zip(doc['products'], doc['reacting_side_chains'], doc['atom_idx']):
#             for filt in deepcopy(regio_filter):
#                 if side_chain == filt['side_chain']:
#                     for tup in filt['threshold_1']:
#                         if atom_idx == tup[0]:
#                             filtered_cands[side_chain]['threshold_1'].append((product, tup[1]))
#                     for tup in filt['threshold_2']:
#                         if atom_idx == tup[0]:
#                             filtered_cands[side_chain]['threshold_2'].append((product, tup[1]))

#         # organize dict
#         for reacting_sc in filtered_cands:
#             for threshold in filtered_cands[reacting_sc]:
#                 prod_heats = list(zip(*sorted(filtered_cands[reacting_sc][threshold], key=lambda x: x[1]))
#                                   ) if filtered_cands[reacting_sc][threshold] else []
#                 filtered_cands[reacting_sc][threshold] = {}
#                 filtered_cands[reacting_sc][threshold]['products'] = prod_heats[0] if prod_heats else []
#                 filtered_cands[reacting_sc][threshold]['heats'] = prod_heats[1] if prod_heats else []

#         # delete entries with no values
#         filtered_cands = clean_results(filtered_cands)

#         results.append((doc['reactant'], filtered_cands))

#     return results

def regiosqm_filter(collection):

    results = []
    for doc in collection:

        # get all unique side_chains that reacted in the reactant and corresponding regioSQM predictions
        unique_sc = list(set([side_chain['side_chain']['ID'] if isinstance(side_chain, dict)
                              else side_chain for side_chain in doc['reacting_side_chains']]))
        unique_sc = [MongoDataBase()['molecules'].find_one({'ID': unique_id})['smiles'] for unique_id in unique_sc]
        # print(unique_sc)
        regio_filter = MongoDataBase()['filters'].find({'smiles': {'$in': unique_sc}})
        # print(regio_filter.count())
        # exit()
        # initialize thresholds for each unique side_chain
        filtered_cands = {'regiosqm': []}

        # place enumerated products in appropriate theshold bin
        for product, side_chain, atom_idx in zip(doc['smiles'], doc['reacting_side_chains'], doc['atom_idxs']):
            for filt in deepcopy(regio_filter):
                if isinstance(side_chain, dict) and side_chain['side_chain']['smiles'] == filt['smiles']:
                    for tup in filt['threshold_1']:
                        if atom_idx == tup[0]:
                            filtered_cands['regiosqm'].append(product)

        # delete entries with no values
        # filtered_cands = clean_results(filtered_cands)

        doc.update(filtered_cands)
        results.append(doc)

    return results


def indole_hetero_atom_filter(collection):
    """
    For each reactant in candidates, apply a filter that selects the enumerated candidates that result from an indole
    nitrogen reacting to close the macrocycle. This filter is generalized to any indole derivative. The filtered
    candidates are stored in a dictionary with the following hierarchy:
        - reacting side chain
            - candidates

    Args:
        collection (iterable): An iterable object containing the candidate data, such as a pymongo cursor on the
            candidates collection

    Returns:
        list: A list containing tuples of the reactant and the sorted dictionary of filtered candidates
    """

    results = []
    for doc in collection:

        # initialize dict based on unique side chains that reacted to close macrocycle
        filtered_cands = {}
        for reacting_sc in list(set(doc['reacting_side_chains'])):
            filtered_cands[reacting_sc] = []

        # apply filter
        for product, side_chain in zip(doc['products'], doc['reacting_side_chains']):
            product = Chem.MolFromSmiles(product)
            if product.HasSubstructMatch(Chem.MolFromSmarts(INDOLE_SMARTS)):
                filtered_cands[side_chain].append(Chem.MolToSmiles(product))

        # delete entries with no values
        filtered_cands = clean_results(filtered_cands)

        results.append((doc['reactant'], filtered_cands))

    return results


def clean_results(filtered_cands):
    """
    Remove any empty fields within the filtered_cands dictionary recursively.

    Args:
        filtered_cands (dict): The dictionary containing the filtered candidates

    Returns:
        dict: The cleaned dictionary
    """

    for key in list(filtered_cands):

        # delete all empty items
        if not filtered_cands[key]:
            del filtered_cands[key]

        # if value is another dictionary, recurse
        elif isinstance(filtered_cands[key], dict):
            filtered_cands[key] = clean_results(filtered_cands[key])

            # delete item if it is empty after cleaning
            if not filtered_cands[key]:
                del filtered_cands[key]

    return filtered_cands


def main():

    parser = argparse.ArgumentParser(description='Apply a set of filters to all candidates within the database and '
                                     'store the filtered candidates in the filtered_molecules collection')
    parser.add_argument('filter', choices=[REGIOSQM_FILTER, INDOLE_FILTER], nargs='+',
                        help='the filter to apply to candidate molecules.')
    parser.add_argument('--db', dest='database', default='molecules',
                        help='the mongoDB database containing the candidates collection')
    parser.add_argument('--host', dest='host', default='localhost',
                        help='the host MongoDB server to connect to')
    parser.add_argument('--port', dest='port', type=int, default=27017,
                        help='the port on host server to connect to')
    parser.add_argument('--show_progress', dest='progress', action='store_false',
                        help='show progress bar (defaults to False)')

    args = parser.parse_args()

    # get candidates
    # db_mol = Database(host=args.host, port=args.port, db=args.database)
    # candidates = db_mol.find_all('candidates', {'_id': 0})

    candidates = MongoDataBase()['molecules'].find({'type': 'candidate'})
    # apply appropriate filters
    if REGIOSQM_FILTER in args.filter:
        results = regiosqm_filter(candidates)
        # for reactant, filtered_cands in tqdm(regiosqm_filter(deepcopy(candidates)), desc='regioSQM: ',
        #                                      disable=args.progress):
        print(len(results))
        pprint(results[:3])
        # if filtered_cands:
        #     Database(db='molecules').insert_filtered_candidates(reactant, REGIOSQM_FILTER, filtered_cands)

    # if INDOLE_FILTER in args.filter:
    #     for reactant, filtered_cands in tqdm(indole_hetero_atom_filter(deepcopy(candidates)), desc='indole: ',
    #                                          disable=args.progress):
    #         if filtered_cands:
    #             Database(db='molecules').insert_filtered_candidates(reactant, INDOLE_FILTER, filtered_cands)


if __name__ == '__main__':
    main()
