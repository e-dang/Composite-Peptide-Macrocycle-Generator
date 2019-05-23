from utils import Database
from copy import deepcopy
from pprint import pprint


def main():
    db_mol = Database(db='molecules')
    db_rxn = Database(db='rxn_templates')

    cands = db_mol.find_all('candidates')
    regio_filter = db_rxn.find_all('regioSQM')

    results = {}
    for candidate in cands:
        print()
        print()
        # print('reactant', candidate['reactant'])
        unique_sc = list(set(candidate['reacting_side_chains']))
        regio_filter = db_rxn.find('regioSQM', {'side_chain': {'$in': unique_sc}}, {'_id': 0})
        # print(type(regio_filter[0]))
        # for x in regio_filter:
        #     print(x)
        # exit()
        rank = {}
        for reacting_sc in unique_sc:
            rank[reacting_sc] = {}
            rank[reacting_sc]['threshold_1'] = []
            rank[reacting_sc]['threshold_2'] = []
        for product, side_chain, atom_idx in zip(candidate['products'], candidate['reacting_side_chains'], candidate['atom_idx']):
            for filt in deepcopy(regio_filter):
                if side_chain == filt['side_chain']:
                    for tup in filt['threshold_1']:
                        if atom_idx == tup[0]:
                            rank[side_chain]['threshold_1'].append((product, tup[1]))
                    for tup in filt['threshold_2']:
                        if atom_idx == tup[0]:
                            rank[side_chain]['threshold_2'].append((product, tup[1]))

        for top_key in rank.keys():
            for inner_key in rank[top_key].keys():
                # rank[top_key][inner_key] = sorted(rank[top_key][inner_key], key=lambda x: x[1])
                # print('rank', rank[top_key][inner_key][0])
                # exit()
                prod_heats = list(zip(*sorted(rank[top_key][inner_key], key=lambda x: x[1]))
                                  ) if rank[top_key][inner_key] else []
                # print(prod_heats)
                rank[top_key][inner_key] = {}
                rank[top_key][inner_key]['products'] = prod_heats[0] if prod_heats else []
                rank[top_key][inner_key]['heats'] = prod_heats[1] if prod_heats else []

        # exit()
        # print(rank)
        # products, heats = zip(*[i for i in sorted(rank.items(), key=lambda x: x[1]) if i[1]])
        # print(items)
        # print(products)
        # print(heats)
        db_mol.insert_regio_filtered_candidates(candidate['reactant'], rank)
        # results[candidate['reactant']] = [i for i in sorted(rank.items(), key=lambda x: x[1]) if i[1]]

    # # pprint(results)
    # with open('/Users/ericdang/Desktop/test2.txt', 'w') as f:
    #     pprint(results, f)


if __name__ == '__main__':
    main()
