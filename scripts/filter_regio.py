import argparse
from copy import deepcopy

from tqdm import tqdm

from utils import Database

REGIOSQM_FILTER = 'regioSQM'


def regioSQM_filter(candidates):
    """
    For each reactant in candidates, the regioSQM filter is applied to the enumerated products resulting from each
    unique side chain reacting to close the macrocycle. The enumerated products are stored in a dictionary with the
    following hierarchy:
        - reactant
        - unique side chain
            - threshold
                - products
                - heats

    Args:
        candidates (iterable): An iterable object containing the candidate data, such as a pymongo cursor on the
            candidates collection.

    Returns:
        list: A list containing tuples of the reactant and sorted dictionary of filtered products and their
            corresponding heats.
    """

    db_rxn = Database(db='rxn_templates')

    results = []
    for candidate in candidates:

        # get all unique side_chains that reacted in the candidate
        unique_sc = list(set(candidate['reacting_side_chains']))
        regio_filter = db_rxn.find('regioSQM', {'side_chain': {'$in': unique_sc}}, {'_id': 0})

        # initialize thresholds for each unique side_chain
        rank = {}
        for reacting_sc in unique_sc:
            rank[reacting_sc] = {}
            rank[reacting_sc]['threshold_1'] = []
            rank[reacting_sc]['threshold_2'] = []

        # place enumerated products in appropriate theshold bin
        for product, side_chain, atom_idx in zip(candidate['products'], candidate['reacting_side_chains'], candidate['atom_idx']):
            for filt in deepcopy(regio_filter):
                if side_chain == filt['side_chain']:
                    for tup in filt['threshold_1']:
                        if atom_idx == tup[0]:
                            rank[side_chain]['threshold_1'].append((product, tup[1]))
                    for tup in filt['threshold_2']:
                        if atom_idx == tup[0]:
                            rank[side_chain]['threshold_2'].append((product, tup[1]))

        # organize dict
        for top_key in rank.keys():
            for inner_key in rank[top_key].keys():
                prod_heats = list(zip(*sorted(rank[top_key][inner_key], key=lambda x: x[1]))
                                  ) if rank[top_key][inner_key] else []
                rank[top_key][inner_key] = {}
                rank[top_key][inner_key]['products'] = prod_heats[0] if prod_heats else []
                rank[top_key][inner_key]['heats'] = prod_heats[1] if prod_heats else []

        results.append((candidate['reactant'], rank))

    return results


def main():

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('filter', choices=[REGIOSQM_FILTER], help='The filter to apply to candidate molecules.')
    parser.add_argument('--db', dest='database', default='molecules',
                        help='The mongoDB database containing candidates.')
    parser.add_argument('--host', dest='host', default='localhost',
                        help='The host MongoDB server to connect to.')
    parser.add_argument('--port', dest='port', type=int, default=27017,
                        help='The port on host server to connect to.')
    parser.add_argument('--show_progress', dest='progress', action='store_false',
                        help='Show progress bar. Defaults to False.')

    args = parser.parse_args()

    db_mol = Database(host=args.host, port=args.port, db=args.database)

    if args.filter == REGIOSQM_FILTER:
        for reactant, rank in tqdm(regioSQM_filter(db_mol.find_all('candidates')), disable=args.progress):
            db_mol.insert_regio_filtered_candidates(reactant, rank)


if __name__ == '__main__':
    main()
