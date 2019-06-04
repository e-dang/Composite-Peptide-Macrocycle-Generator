from rdkit import Chem
from rdkit.Chem import AllChem
from pathlib import Path
import argparse
from utils import Database
from tqdm import tqdm
from rdkit.Chem import rdmolfiles
from time import time

COLLECTION = 'filtered_candidates'
FILTER_REGIO = 'regioSQM'
FILTER_INDOLE = 'indole'


def generate_conformers(filtered_candidates, num_confs):

    for candidate in tqdm(filtered_candidates, desc='Candidates: ', total=filtered_candidates.count()):
        products = get_products(candidate)
        yield embed_conformers(products, num_confs)


def get_products(candidate):

    if candidate['filter'] == FILTER_REGIO:
        return get_regiosqm_products(candidate)
    if candidate['filter'] == FILTER_INDOLE:
        return get_indole_products(candidate)


def get_regiosqm_products(candidate):

    attr = 'reacting_side_chains'
    products = []
    for side_chain in candidate[attr]:
        for threshold in candidate[attr][side_chain]:
            for product in candidate[attr][side_chain][threshold]['products']:
                products.append(product)

    return products


def get_indole_products(candidate):

    attr = 'reacting_side_chains'
    products = []
    for side_chain in candidate[attr]:
        for product in candidate[attr][side_chain]:
            products.append(product)

    return products


def embed_conformers(products, num_confs):

    # convert to list if not already a list
    if not isinstance(products, list):
        products = [products]

    # convert all SMILES strings to mol and add hydrogens
    products = [Chem.MolFromSmiles(mol) for mol in products if isinstance(mol, str)]
    [AllChem.AddHs(mol) for mol in products]

    # embed molecule
    for mol in products:
        conf_ids = AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, numThreads=0, randomSeed=int(time()))
        return (mol, conf_ids)


def optimize_conformers(conformers):

    for mol, conf_ids in conformers:
        convergence, energy = zip(*AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=1000, numThreads=0))

        # sort confomers based on energy for proper alignment
        sorted_confs = sorted(zip(energy, conf_ids), key=lambda x: x[0])
        energy, conf_ids = zip(*sorted_confs)

        # align
        rms = []
        AllChem.AlignMolConformers(mol, confIds=conf_ids, maxIters=1000, RMSlist=rms)

        # check if every conformer converged
        if len(set(convergence)) != 1:
            analyze_convergence(mol, conf_ids, convergence)

        yield (mol, conf_ids, rms, convergence, energy)


def analyze_convergence(mol, conf_ids, convergence):

    for indicator, conf_id in zip(convergence, conf_ids):
        if indicator != 0:
            print('Did not converge!')
            print('Candidate:', Chem.MolToSmiles(mol))
            print('Conformer ID:', conf_id)


def save_conformers(confs, **kwargs):

    db = Database(host=kwargs['host'], port=kwargs['port'], db=kwargs['database'])
    for mol, conf_ids, rms, convergences, energies in confs:
        db.insert_conformers(Chem.MolToSmiles(mol), mol.ToBinary(), convergences, energies)

        if Chem.MolToSmiles(mol) == 'O=C1CCc2cccc(c2)/C=C/Cc2nn3c(c2CCC[C@H](C(=O)N[C@@H](CCCc2csc4ncsc24)C(=O)NC[C@@H](CCc2csc4ncsc24)C(=O)O)CN1)NC=N3':
            rdmolfiles.MolToPDBFile(mol, 'GT.pdb')
            print(rms)
            print(energies)
            # for i, conf in enumerate(Chem.rdchem.Mol.GetConformers(mol)):
            #     print(rms)
            #     rdmolfiles.MolToPDBFile(mol, f'GT{i}.pdb', confId=i)
            exit()


def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--num_confs', default=1000, type=int,
                        help='The number of conformers to generate for each filtered candidate')
    parser.add_argument('--optimize', dest='optimize', action='store_false',
                        help='Toggles whether to run MMFF94 optimization on the conformers')
    parser.add_argument('--db', dest='database', default='molecules',
                        help='The mongoDB database to connect to')
    parser.add_argument('--host', dest='host', default='localhost',
                        help='The host MongoDB server to connect to')
    parser.add_argument('--port', dest='port', type=int, default=27017,
                        help='The port on host server to connect to')

    args = parser.parse_args()

    db_mol = Database(host=args.host, port=args.port, db=args.database)
    filtered_candidates = db_mol.find_all(COLLECTION)

    confs = generate_conformers(filtered_candidates, args.num_confs)

    if args.optimize:
        opt_confs = optimize_conformers(confs)
        save_conformers(opt_confs, **{'host': args.host, 'port': args.port, 'database': args.database})

    return True


if __name__ == '__main__':
    main()
