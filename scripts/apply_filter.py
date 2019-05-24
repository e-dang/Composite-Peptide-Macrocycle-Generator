import argparse
from copy import deepcopy

from rdkit import Chem
from tqdm import tqdm

from utils import Database

REGIOSQM_FILTER = 'regioSQM'
INDOLE_FILTER = 'indole'

INDOLE = '*/C=C/Cn1aaa2aaaaa21'


def regioSQM_filter(collection):
    """
    For each doc in the collection, the regioSQM filter is applied to the enumerated candidates resulting from each
    unique side chain reacting to close the macrocycle. The filtered candidates are stored in a dictionary with the
    following hierarchy:
        - reacting side chain
            - threshold
                - candidates
                - heats

    Args:
        collection (iterable): An iterable object containing the candidate data, such as a pymongo cursor on the
            candidates collection

    Returns:
        list: A list containing tuples of the reactant and the sorted dictionary of filtered candidates
    """

    db_rxn = Database(db='rxn_templates')

    results = []
    for doc in collection:

        # get all unique side_chains that reacted in the reactant and corresponding regioSQM predictions
        unique_sc = list(set(doc['reacting_side_chains']))
        regio_filter = db_rxn.find('regioSQM', {'side_chain': {'$in': unique_sc}}, {'_id': 0})

        # initialize thresholds for each unique side_chain
        filtered_prods = {}
        for reacting_sc in unique_sc:
            filtered_prods[reacting_sc] = {}
            filtered_prods[reacting_sc]['threshold_1'] = []
            filtered_prods[reacting_sc]['threshold_2'] = []

        # place enumerated products in appropriate theshold bin
        for product, side_chain, atom_idx in zip(doc['products'], doc['reacting_side_chains'], doc['atom_idx']):
            for filt in deepcopy(regio_filter):
                if side_chain == filt['side_chain']:
                    for tup in filt['threshold_1']:
                        if atom_idx == tup[0]:
                            filtered_prods[side_chain]['threshold_1'].append((product, tup[1]))
                    for tup in filt['threshold_2']:
                        if atom_idx == tup[0]:
                            filtered_prods[side_chain]['threshold_2'].append((product, tup[1]))

        # organize dict
        for top_key in filtered_prods:
            for inner_key in filtered_prods[top_key].keys():
                prod_heats = list(zip(*sorted(filtered_prods[top_key][inner_key], key=lambda x: x[1]))
                                  ) if filtered_prods[top_key][inner_key] else []
                filtered_prods[top_key][inner_key] = {}
                filtered_prods[top_key][inner_key]['products'] = prod_heats[0] if prod_heats else []
                filtered_prods[top_key][inner_key]['heats'] = prod_heats[1] if prod_heats else []

        results.append((doc['reactant'], filtered_prods))

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

        # get all unique side_chains that reacted in the reactant
        unique_sc = list(set(doc['reacting_side_chains']))

        # initialize dict
        filtered_prods = {}
        for reacting_sc in unique_sc:
            filtered_prods[reacting_sc] = []

        # apply filter
        for product, side_chain in zip(doc['products'], doc['reacting_side_chains']):
            product = Chem.MolFromSmiles(product)
            if product.HasSubstructMatch(Chem.MolFromSmarts(INDOLE)):
                filtered_prods[side_chain].append(Chem.MolToSmiles(product))

        results.append((doc['reactant'], filtered_prods))

    return results


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
    db_mol = Database(host=args.host, port=args.port, db=args.database)
    candidates = db_mol.find_all('candidates', {'_id': 0})

    # apply appropriate filters
    if REGIOSQM_FILTER in args.filter:
        for reactant, filtered_prods in tqdm(regioSQM_filter(candidates), desc='regioSQM: ', disable=args.progress):
            Database(db='molecules').insert_filtered_candidates(reactant, REGIOSQM_FILTER, filtered_prods)

    if INDOLE_FILTER in args.filter:
        for reactant, filtered_prods in tqdm(indole_hetero_atom_filter(candidates), desc='indole: ',
                                             disable=args.progress):
            Database(db='molecules').insert_filtered_candidates(reactant, INDOLE_FILTER, filtered_prods)


if __name__ == '__main__':
    main()
