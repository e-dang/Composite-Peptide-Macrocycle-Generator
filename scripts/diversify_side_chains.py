import argparse
import json
from copy import copy
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

from utils import read_mols

SIDE_CHAIN_MAP_NUM = 1
CONNECTION_MAP_NUM = 2


def alternate_connection_point_test(mol, connection):

    mols = set()
    attached = set()
    connection = Chem.MolFromSmarts(connection)

    # check if connecting atom is atom mapped
    map_nums = [atom.GetAtomMapNum() for atom in connection.GetAtoms()]
    if CONNECTION_MAP_NUM not in map_nums:
        print('Need to specifiy connecting atom with atom map number')
        raise SystemExit

    # try to make attachment at each atom
    for atom in mol.GetAtoms():

        # detetmine atom eligibility
        atom_idx = None
        found = False
        if atom.GetSymbol() == 'C' and atom.GetTotalNumHs() != 0 and atom.GetTotalNumHs() < 3:
            atom.SetAtomMapNum(SIDE_CHAIN_MAP_NUM)
            atom_idx = atom.GetIdx()
            attached.add(atom_idx)
            found = True
        elif (atom.GetSymbol() == 'N' or atom.GetSymbol() == 'O' or atom.GetSymbol() == 'S') and atom.GetTotalNumHs() != 0:
            atom.SetAtomMapNum(SIDE_CHAIN_MAP_NUM)
            atom_idx = atom.GetIdx()
            attached.add(atom_idx)
            found = True
        elif atom.GetIdx() in attached or not found:
            continue

        # prepare for attachment
        combo = Chem.RWMol(Chem.CombineMols(mol, connection))

        # find reacting atoms on combo mol and reset atom map numbers
        mol_atom = None
        conn_atom = None
        for combo_atom in combo.GetAtoms():
            if combo_atom.GetAtomMapNum() == SIDE_CHAIN_MAP_NUM:
                mol_atom = combo_atom.GetIdx()
                combo_atom.SetAtomMapNum(0)
                Chem.Mol.GetAtomWithIdx(mol, atom_idx).SetAtomMapNum(0)
            elif combo_atom.GetAtomMapNum() == CONNECTION_MAP_NUM:
                conn_atom = combo_atom.GetIdx()
                combo_atom.SetAtomMapNum(0)

        # create bond
        combo.AddBond(mol_atom, conn_atom, order=Chem.rdchem.BondType.SINGLE)

        # fix hydrogen counts
        atom_react = combo.GetAtomWithIdx(mol_atom)
        if atom_react.GetSymbol() == 'N' or atom_react.GetSymbol() == 'O' or atom_react.GetSymbol() == 'S':
            atom_react.SetNumExplicitHs(0)
        elif atom_react.GetSymbol() == 'C' and Chem.Atom.GetTotalNumHs(atom_react) > 0:
            atom_react.SetNumExplicitHs(Chem.Atom.GetTotalNumHs(atom_react) - 1)

        try:
            Chem.SanitizeMol(combo)
            mols.add(Chem.MolToSmiles(combo))
        except ValueError:
            print("Can't sanitize mol:", Chem.MolToSmiles(mol))

    return mols


def transform_connection_point(mols, new_connection, conn_idx, old_connection='[CH3]'):

    new_smiles = set()
    for mol in mols:
        new_mols = Chem.ReplaceSubstructs(Chem.MolFromSmiles(mol), Chem.MolFromSmarts(old_connection),
                                          Chem.MolFromSmarts(new_connection), replacementConnectionPoint=conn_idx)
        for new_mol in new_mols:
            Chem.SanitizeMol(new_mol)
            new_smiles.add(Chem.MolToSmiles(new_mol))

    return new_smiles


def accumulate_mols(mols, collection, parent, modifications, group):

    for smiles in mols:
        doc = {}
        doc['side_chain'] = smiles
        doc['parent'] = parent
        doc['modification'] = modifications
        doc['group'] = group
        collection.append(doc)


def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i', '--in', dest='input', nargs='+',
                        default=['side_chains_likely1.sdf'], help='The sdf file containing monomer side chains.')
    parser.add_argument('-o', '--out', dest='out', default='side_chains.json', help='The output json file.')
    parser.add_argument('-fi', '--fin', dest='fp_in', default='chemdraw/pre_monomer/',
                        help='The input filepath relative to script')
    parser.add_argument('-fo', '--fout', dest='fp_out', default='smiles/pre_monomer/',
                        help='The ouput filepath relative to script')

    args = parser.parse_args()

    # get side_chain group name from input file name
    groups = [name.split('_')[-1].split('.')[0] for name in args.input]

    # get absolute filepath to input and output files
    base_path = Path(__file__).resolve().parents[1]
    fp_in = [str(base_path / args.fp_in / file) for file in args.input]
    fp_out = str(base_path / args.fp_out / args.out)

    # # get mols from all files
    # mols = read_multiple_sdf(fp_in)

    # diversify
    collection = []
    for fp, group in zip(fp_in, groups):
        for mol in read_mols(fp):
            unique_methyl_mols = alternate_connection_point_test(mol, '[CH3:2]')
            unique_ethyl_mols = alternate_connection_point_test(mol, '[CH3][CH2:2]')
            unique_propyl_mols = alternate_connection_point_test(mol, '[CH3][CH2][CH2:2]')
            accumulate_mols(unique_methyl_mols, collection, Chem.MolToSmiles(mol), [0, 3], group)
            accumulate_mols(unique_ethyl_mols, collection, Chem.MolToSmiles(mol), [1, 3], group)
            accumulate_mols(unique_propyl_mols, collection, Chem.MolToSmiles(mol), [2, 3], group)

    # write data to json file
    with open(fp_out, 'w') as f:
        json.dump(collection, f)


if __name__ == '__main__':
    main()
