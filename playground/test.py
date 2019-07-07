from copy import deepcopy
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import rdChemReactions
from rdkit.Chem import rdmolfiles
from macrocycles.utils import write_mols, read_mols, create_monomer_requirements, kekulize_smiles, conformers_to_pdb, PROJECT_DIR, MongoDataBase, DataInitializer, MongoParams, Base
from macrocycles.config import COLLECTIONS, DATA_DIR, COL1, SCM_DOC_TYPE
import json
from itertools import islice
from pymongo import MongoClient, ASCENDING
from collections import OrderedDict
from bson import ObjectId, json_util
import os
import numpy as np
from copy import deepcopy
from macrocycles.generate_peptides import PeptideGenerator
from macrocycles.merge_template_peptide import TPHybridGenerator
# rxn = AllChem.ReactionFromSmarts(
#     '([*:3]Cc1cccc(/C=C/[C:1]OC(=O)OC(C)(C)C)c1.c1nn2n[cH:2]c(C[*:4])c2[nH]1)>>C(=C/[C:1][c:2]1nn2nc[nH]c2c1C[*:4])\\c1cccc(C[*:3])c1')
# mol = Chem.MolFromSmarts('CC(C)(C)OC(=O)O[CH2:1]/C=C/C1=CC(C[*:3])=CC=C1')
# mol = Chem.MolFromSmiles(
#     'CC(C)(C)OC(=O)OC/C=C/C1=CC=CC(CCC(=O)NC[C@H](CCCC2=C3[NH]C=NN3N=C2)C(=O)N[C@@H](CCCC2=CSC3=C2SC=N3)C(=O)NC[C@@H](CCC2=CSC3=C2SC=N3)C(=O)O)=C1')
# db = Database(db='molecules')
# cur = db.find('candidates', {
#               'reactant': 'CC(C)(C)OC(=O)OC/C=C/c1cccc(CCC(=O)NC[C@H](CCCc2cnn3nc[nH]c23)C(=O)N[C@@H](CCCc2csc3ncsc23)C(=O)NC[C@@H](CCc2csc3ncsc23)C(=O)O)c1'})
# for doc in cur:
#     for mol in doc['products']:
# mol = Chem.MolFromSmiles(
#     'N[C@@H](Cc1csc2ncsc12)C(=O)O')
# Chem.Kekulize(mol)
# print(Chem.MolToSmiles(mol, kekuleSmiles=True))
# print(Chem.MolToSmiles(mol, kekuleSmiles=True))
# Draw.MolToImage(mol, size=(1000, 1000), includeAtomNumbers=True).show()
# Draw.MolToImage(mol, size=(500, 500)).show()
# Chem.Kekulize(mol)
# print(Chem.MolToSmiles(mol, kekuleSmiles=True))
# Draw.MolToImage(mol, size=(500, 500)).show()
# Draw.ReactionToImage(rxn, subImgSize=(500, 500)).show()

# products = rxn.RunReactants((reactant,))
# print(len(products))
# prods = set()
# for product in products:
#     Chem.SanitizeMol(product[0])
#     prods.add(Chem.MolToSmiles(product[0]))
# for p in prods:
#     Draw.MolToImage(Chem.MolFromSmiles(p), size=(500, 500)).show()


# mol = Chem.MolFromSmiles(
#     'O=C1CCc2cc(/C=C/Cc(cc3)cc4c3[nH]cc4C[C@H](N1)C(N[C@H](C(N[C@H](C(O)=O)Cc5ccc(O)cc5)=O)Cc(c[nH]6)c7c6cccc7)=O)ccc2')

# major_prod = 'O=C1CCc2cc(/C=C/Cc(cc3)cc4c3[nH]cc4C[C@H](NC([C@@H](N1)Cc(c[nH]5)c6c5cccc6)=O)C(N[C@H](C(O)=O)Cc7ccc(O)cc7)=O)ccc2'
# major_prod2 = 'O=C1CCc2cc(/C=C/Cc(cc3)cc4c3[nH]cc4C[C@H](N1)C(N[C@H](C(N[C@H](C(O)=O)Cc5ccc(O)cc5)=O)Cc(c[nH]6)c7c6cccc7)=O)ccc2'
# linear = 'CC(C)(OC(OC/C=C/c1cc(CCC(N[C@H](C(N[C@H](C(N[C@H](C(O)=O)Cc2ccc(O)cc2)=O)Cc(c[nH]3)c4c3cccc4)=O)Cc(c[nH]5)c6c5cccc6)=O)ccc1)=O)C'
# Chem.Kekulize(mol)
# print(Chem.MolToSmiles(mol, kekuleSmiles=True))
# print(Chem.MolToSmiles(mol))
# tmp = AllChem.Compute2DCoords(mol)
# AllChem.COmp
# print(tmp)
# AllChem.Ge
# AllChem.AddHs(mol)
# AllChem.EmbedMolecule(mol)
# confs = AllChem.EmbedMultipleConfs(mol, numConfs=10)
# converge, energy = zip(*AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=200))
# bi = mol.ToBinary()
# m = Chem.Mol(bi)
# print(m.GetNumConformers())
# for conf, eng in zip(confs, energy):
# mol.SetProp('_Name', 'Conf_ID:' + str(conf) + '\nEnergy:' + str(eng))
# block = Chem.MolToMolBlock(mol, confId=conf)
# print([block])
# mol = Chem.MolFromMolBlock(block)
# exit()
#     break
#     print('\n\n')
# print(AllChem.MMFFOptimizeMolecule(mol, maxIters=1000))
# ff = AllChem.MMFFGetMoleculeForceField(mol)
# print(ff.CalcEnergy())
# ls = []
# AllChem.AlignMolConformers(mol, RMSlist=ls)
# print(ls)
# print(len(confs))
# rdmolfiles.MolToPDBFile(mol, 'WWY_major_prod2_optimized.pdb')
# print('done')
# for i, molec in enumerate(confs):

# rdmolfiles.MolToPDBFile(mol, 'test' + str(i) + '.pdb')
# print(AllChem.MMFFOptimizeMolecule(mol, maxIters=1000))
# rdmolfiles.MolToPDBFile(mol, 'test.pdb')
# writer = Chem.SDWriter('test.sdf')
# writer.write(mol)
# Draw.MolToImage(mol).show()

# with open('/Users/ericdang/Documents/UCLA_Research/macrocycles/smiles/monomers/modified_prolines.json', 'w') as fout:
#     mols = read_mols('/Users/ericdang/Documents/UCLA_Research/macrocycles/chemdraw/monomers/modified_prolines.sdf')
#     collection = []
#     for mol in mols:
#         doc = {}
#         doc['monomer'] = Chem.MolToSmiles(mol)
#         doc['type'] = 'alpha_amino_acid'
#         doc['side_chain'] = ''
#         collection.append(doc)
#     json.dump(collection, fout)

# with open('/Users/ericdang/Documents/UCLA_Research/macrocycles/smiles/monomers/modified_prolines_CCW.json', 'r') as f:
#     for doc in json.load(f):
#         mol = Chem.MolFromSmiles(doc['monomer'])
#         print(Chem.MolToSmiles(mol))
# patt = Chem.MolFromSmarts('NCC(=O)O')
# matches = mol.GetSubstructMatches(patt, useChirality=False)
# for pair in matches:
#     for atom_idx in pair:
#         atom = Chem.Mol.GetAtomWithIdx(mol, atom_idx)
#         if atom.GetSymbol() == 'N' and Chem.Atom.GetTotalNumHs(atom) != 0:
#             atom.SetAtomMapNum(1)
# print(Chem.MolToSmiles(mol))
# Draw.MolToImage(mol).show()
# x = input()


# mol = Chem.RWMol(Chem.MolFromSmiles('O=C1O[C@@H]2CN[C@H]1C2'))
# patt = Chem.MolFromSmarts('NCC(=O)O')
# match = mol.GetSubstructMatches(patt)
# for pair in match:
#     for atom_idx in pair:
#         atom = Chem.Mol.GetAtomWithIdx(mol, atom_idx)
#         # print(Chem.Atom.GetHybridization(atom))
#         # if Chem.Atom.GetHybridization(atom) == Chem.rdchem.HybridizationType.SP2:
#         #     print(atom.GetSymbol())
#         if atom.GetSymbol() == 'O' and Chem.Atom.GetTotalNumHs(atom) == 1:
#             mol.RemoveAtom(atom_idx)
#         elif atom.GetSymbol() == 'C' and Chem.Atom.GetTotalNumHs(atom) == 0 and Chem.Atom.GetHybridization(atom) == Chem.rdchem.HybridizationType.SP2:
#             atom.SetAtomMapNum(2)

# Draw.MolToImage(mol).show()
# print(Chem.MolToSmiles(mol))


# mol = Chem.MolFromSmiles(
#     'O=C1CCc2cccc(c2)/C=C/Cc2cc(ccc2O)C[C@@H](C(=O)O)NC(=O)[C@H](Cc2c[nH]c3ccccc23)NC(=O)[C@H](Cc2c[nH]c3ccccc23)N1')
# Chem.Kekulize(mol)
# print(Chem.MolToSmiles(mol, kekuleSmiles=True))

# with open('smiles.txt', 'r') as f:
#     smiles = [line.strip() for line in f.readlines()]
#     kekulize_smiles(smiles)

# rxn = AllChem.ReactionFromSmarts(
#     '([*:3]Cc1cccc(/C=C/[C:1]OC(=O)OC(C)(C)C)c1.[*:5][C@@H]([NH:7][*:4])C[c:2]1[c:6][nH]c2ccccc12)>>C(=C/[C:1][C:2]12C[C@@H]([*:5])[N:7]([*:4])[C:6]1Nc1ccccc12)\c1cccc(C[*:3])c1')
# rxn = AllChem.ReactionFromSmarts(
#     '([*:3]Cc1cccc(/C=C/[C:1]OC(=O)OC(C)(C)C)c1.[*:5][C@@H]([NH:7][*:4])C[c:2]1[c:6][nH]c2ccccc12)>>C(=C/[C:1][C:2]12C[C@@H]([*:5])[N:7]([*:4])[C:6]1Nc1ccccc12)\\c1cccc(C[*:3])c1')
# mol = Chem.MolFromSmiles(
#     "CC(C)(C)OC(=O)OC/C=C/c1cccc(CCC(=O)N[C@@H](Cc2c[nH]c3ccccc23)C(=O)N[C@@H](Cc2c[nH]c3ccccc23)C(=O)N[C@@H](Cc2ccc(O)cc2)C(=O)O)c1")
# Chem.Kekulize(mol)
# print(Chem.MolToSmiles(mol, kekuleSmiles=True))
# Draw.MolToImage(mol, size=(500, 500)).show()
# Draw.ReactionToImage(rxn, subImgSize=(500, 500)).show()
# prods = rxn.RunReactants((mol,))
# print(len(prods))
# prod = set()
# for m in prods:
#     Chem.SanitizeMol(m[0])
#     prod.add(Chem.MolToSmiles(m[0]))

# for m in prod:
#     print(m)
#     Draw.MolToImage(Chem.MolFromSmiles(m, ), size=(500, 500)).show()

# conformers_to_pdb(
#     'O=C1CCc2cccc(c2)/C=C/Cn2cc(c3ccccc32)C[C@@H](C(=O)N[C@@H](Cc2c[nH]c3ccccc23)C(=O)N[C@@H](Cc2ccc(O)cc2)C(=O)O)N1')
# print(Chem.rdMolDescriptors.CalcNumRotatableBonds(Chem.MolFromSmiles(
#     'O=C1CCc2cccc(c2)/C=C/Cc2[nH]c3ccccc3c2C[C@@H](C(=O)N[C@@H](Cc2ccc(O)cc2)C(=O)O)NC(=O)[C@H](Cc2c[nH]c3ccccc23)N1'), True))
# print(Chem.MolToSmiles(Chem.MolFromSmiles('C1=CC=CN1')))

# mol = Chem.MolFromSmiles(
#     'C#CCCC[C@@](CC(=O)N[C@H](CCCc1cccc2c(=O)ncsc12)CC(=O)N[C@@H](CCc1cc2[nH]ccc2cc1F)CC(=O)N[C@@H](CCC(N)=O)C(=O)O)(Cc1cccc(/C=C/COC(=O)OCC(C)C)c1)C(OCC)OCC')
# Chem.Kekulize(mol)
# print(Chem.MolToSmiles(mol, kekuleSmiles=True))
# from time import time, sleep


mol1 = Chem.MolFromSmiles('O=C1N2C(N=CC=C2)=C(C)C=C1')
mol2 = Chem.MolFromSmiles('Cc1ccccc1O')
binstr1 = mol1.ToBinary()
# Chem.Kekulize(mol)
binstr2 = mol2.ToBinary()
print(binstr1 == binstr2)
print(type(binstr1))
print(Chem.MolToSmiles(mol1))
print(Chem.MolToSmiles(mol2))
# print(binstr1)
# print(binstr2)
exit()


generator = PeptideGenerator(f_in='', f_out='')
generator.generate_peptide_from_ids(['a1AA5', 'a1AA5', 'a1DG3'])
print(generator.result_data)
db = MongoDataBase()
db.insert('molecules', generator.result_data)
exit()

# s = set()
# rxn = AllChem.ReactionFromSmarts(
#     '(O=C(C(C[CH:3]=O)[*:2])[NH:4]C(Cc1c2ccccc2[nH][cH:5]1)[*:1])>>O=C1C([*:2])C[CH:3]2[N:4]1C([*:1])Cc1c3ccccc3[nH][c:5]12')
# Draw.ReactionToImage(rxn, subImgSize=(500, 500)).show()
# mol = Chem.MolFromSmiles(
#     'C#CCCC[C@@](C=O)(CC(=O)N[C@@H](Cc1c[nH]c2ccccc12)C(=O)N[C@@H](Cc1c[nH]c2ccccc12)C(=O)N[C@@H](Cc1ccc(O)cc1)C(=O)O)Cc1cccc(/C=C/COC(=O)OC(C)(C)C)c1')
# products = rxn.RunReactants((mol,))
# print(len(products))
# for prod in products:
#     for x in prod:
#         s.add(Chem.MolToSmiles(x))
# for x in s:
#     Draw.MolToImage(Chem.MolFromSmiles(x), size=(500, 500)).show()
# exit()

db = MongoDataBase()
db.setup(clear=True)
db['reactions'].insert([
    {
        "ID": "1",
        "type": "reaction",
        "smarts": "(O=C(C(C[CH:4]=O)[*:2])[NH:3]C(Cc1cc[nH][cH:5]1)[*:1])>>O=C1C([*:2])C[CH:4]2[N:3]1C([*:1])Cc1cc[nH][c:5]12",
        "side_chain": "Cc1c[nH]cc1",
        "rxn_vector": [
            [
                14,
                8
            ],
            [
                14,
                7
            ],
            [
                15,
                0
            ]
        ]
    },
    {
        "ID": "2",
        "type": "reaction",
        "smarts": "(O=C(C(C[CH:3]=O)[*:2])[NH:4]C(Cc1c2ccccc2[nH][cH:5]1)[*:1])>>O=C1C([*:2])C[CH:3]2[N:4]1C([*:1])Cc1c3ccccc3[nH][c:5]12",
        "side_chain": "AA5",
        "rxn_vector": [
            [
                18,
                12
            ],
            [
                18,
                4
            ],
            [
                19,
                0
            ]
        ]
    },
    {
        "ID": "3",
        "type": "reaction",
        "smarts": "(O=C(C(C[CH:4]=O)[*:2])[NH:3]C(Cc1ccs[cH:5]1)[*:1])>>O=C1C([*:2])C[CH:4]2[N:3]1C([*:1])Cc1ccs[c:5]12",
        "side_chain": "Cc1ccsc1",
        "rxn_vector": [
            [
                14,
                8
            ],
            [
                14,
                7
            ],
            [
                15,
                0
            ]
        ]
    },
    {
        "ID": "pyroloindoline",
        "type": "reaction",
        "smarts": "([*:3]Cc1cccc(/C=C/[C:1]OC(=O)OC(C)(C)C)c1.[*:5][C@@H]([NH:7][*:4])C[c:2]1[c:6][nH]c2ccccc12)>>C(=C/[C:1][C:2]12C[C@@H]([*:5])[N:7]([*:4])[C:6]1Nc1ccccc12)\c1cccc(C[*:3])c1",
        "side_chain": "AA5",
        "rxn_vector": []
    },
    {
        "ID": "4",
        "type": "reaction",
        "smarts": "(O=C(CC([*:2])([*:3])[CH:4]=O)[NH:5]C(Cc1c2ccccc2[nH][cH:6]1)[*:1])>>O=C1CC([*:2])([*:3])[CH:4]2[N:5]1C([*:1])Cc1c3ccccc3[nH][c:6]12",
        "side_chain": "AA5",
        "rxn_vector": [
            [
                19,
                12
            ],
            [
                19,
                4
            ],
            [
                20,
                0
            ]
        ]
    }
])
# exit()
db['reactions'].insert([
    {
        "ID": "t1",
        "type": "template",
        "smarts": "[*:3]Cc1cccc(/C=C/[CH2:1]OC(=O)OC(C)(C)C)c1",
        "matches": [
            "t1",
            "t3"
        ],
        "atom_idx": 9
    },
    {
        "ID": "t2",
        "type": "template",
        "smarts": "[*:3]Cc1cc(/C=C/[CH2:1]OC(=O)OC(C)(C)C)ccc1F",
        "matches": [
            "t2"
        ],
        "atom_idx": 7
    }
])

file_loader = DataInitializer('pre_monomer/side_chains_likely1.sdf', 'pre_monomer/side_chains_test')
file_loader.load_parent_side_chains()
file_loader = DataInitializer('pre_monomer/monomer_backbone.sdf', 'pre_monomer/backbones_test')
file_loader.load_backbones()
file_loader = DataInitializer('templates.sdf', 'peptides/templates')
file_loader.load_templates()
exit()


# s = set()
# rxn = AllChem.ReactionFromSmarts(
#     '(CC(C)(C)OC(=O)O[CH2:1]/C=C/c1cccc(C[*:3])c1.c1oc[cH:2]c1C[*:4])>>C(=C/[CH2:1][c:2]1cocc1C[*:4])\\c1cccc(C[*:3])c1')
# mol = Chem.MolFromSmiles(
#     'CC(C)(C)OC(=O)OC/C=C/c1cccc(CCC(=O)N[C@H](Cc2ccoc2)C(=O)N[C@H](Cc2ccoc2)C(=O)N[C@H](Cc2ccoc2)C(=O)O)c1')
# products = rxn.RunReactants((mol,))
# for prod in products:
#     for x in prod:
#         s.add(Chem.MolToSmiles(x))

# for x in s:
#     Draw.MolToImage(Chem.MolFromSmiles(x), size=(500, 500)).show()
# Draw.ReactionToImage(rxn, subImgSize=(500, 500)).show()

mol = Chem.MolFromSmiles('*C(Cc1c[nH]c2ccccc12)NC(=O)CC(*)(*)C=O')
print(Chem.MolToSmiles(mol))
Draw.MolToImage(mol, size=(500, 500), includeAtomNumbers=True).show()


def run_reaction_vector(mol, vector=None):

    # for atom in mol.GetAtoms():
    #     map_num = atom.GetAtomMapNum()
    #     if map_num != 0:
    mol = mol.split('*')
    for i, substr in enumerate(mol):
        if i == 0:
            mol = substr
        else:
            mol = f'[*:{i}]'.join([mol, substr])

    print(mol)
    mol = Chem.MolFromSmiles(mol)

    idxs = [idx for idx_tup in vector if 0 not in idx_tup for idx in idx_tup]
    idxs = set(idxs)
    print(idxs)
    for i, idx in enumerate(idxs, start=i + 1):
        mol.GetAtomWithIdx(idx).SetAtomMapNum(i)

    old_mol = deepcopy(mol)
    mol = Chem.RWMol(mol)
    idxs = set()
    for idx_tup in vector:
        idx1, idx2 = idx_tup
        idxs.add(idx1)
        idxs.add(idx2)

        if 0 in idx_tup:
            idx = idx_tup[0] if idx_tup[0] != 0 else idx_tup[1]
            mol.RemoveAtom(idx)
        else:
            idx1, idx2 = idx_tup
            idxs.add(idx1)
            idxs.add(idx2)
            mol.AddBond(idx1, idx2, Chem.rdchem.BondType.SINGLE)

    # idxs = [idx for idx_tup in vector for idx in idx_tup if idx != 0]
    # idxs = set(idxs)

    # mapping = {0: 0}
    # for i, idx in enumerate(idxs, start=1):
    #     old_mol.GetAtomWithIdx(idx).SetAtomMapNum(i)
    #     mol.GetAtomWithIdx(idx).SetAtomMapNum(i)
    reaction = '(' + Chem.MolToSmiles(old_mol) + ')>>' + Chem.MolToSmiles(mol)

    # wild_cards = [i for i, char in enumerate(reaction) if char == '*']
    # for i in wild_cards:
    #     reaction
    return reaction


vector = [(19, 12), (19, 4), (20, 0)]
reaction = run_reaction_vector('*C(Cc1c[nH]c2ccccc12)NC(=O)CC(*)(*)C=O', vector=vector)
# reaction = '(' + Chem.MolToSmiles(mol) + ')>>' + Chem.MolToSmiles(product)
print(reaction)
rxn = AllChem.ReactionFromSmarts(
    '(O=C(C(C[CH:3]=O)[*:2])[NH:4]C(Cc1c2ccccc2[nH][cH:5]1)[*:1])>>O=C1C([*:2])C[CH:3]2[N:4]1C([*:1])Cc1c3ccccc3[nH][c:5]12')
Draw.ReactionToImage(rxn, subImgSize=(500, 500)).show()
s = set()
# rxn = AllChem.ReactionFromSmarts(
#     '(CC(C)(C)OC(=O)O[CH2:1]/C=C/c1cccc(C[*:3])c1.c1oc[cH:2]c1C[*:4])>>C(=C/[CH2:1][c:2]1cocc1C[*:4])\\c1cccc(C[*:3])c1')
mol = Chem.MolFromSmiles(
    'CC(C)(C)OC(=O)OC/C=C/c1ccc(F)c(C[C@@H](CC=O)C(=O)N[C@H](Cc2c[nH]c3ccccc23)C(=O)N[C@H](Cc2c[nH]c3ccccc23)C(=O)N[C@H](Cc2ccc(O)cc2)C(=O)O)c1')
# mol = Chem.MolFromSmiles(
#     'C#CCCC[C@@](C=O)(CC(=O)N[C@@H](Cc1c[nH]c2ccccc12)C(=O)N[C@@H](Cc1c[nH]c2ccccc12)C(=O)N[C@@H](Cc1ccc(O)cc1)C(=O)O)Cc1cccc(/C=C/COC(=O)OC(C)(C)C)c1')
Draw.MolToImage(mol, size=(500, 500)).show()
products = rxn.RunReactants((mol,))
print(len(products))
for prod in products:
    for x in prod:
        s.add(Chem.MolToSmiles(x))
for x in s:
    Draw.MolToImage(Chem.MolFromSmiles(x), size=(500, 500)).show()

# monomers = [{
#     'smiles': 'N[C@H](Cc1c[nH]c2ccccc12)C(=O)O',
#     'backbone': 'alpha',
#     'required': True
# }, {
#     'smiles': 'N[C@H](Cc1c[nH]c2ccccc12)C(=O)O',
#     'backbone': 'alpha',
#     'required': True
# }, {'smiles': 'N[C@H](Cc1ccc(O)cc1)C(=O)O',
#     'backbone': 'alpha',
#     'required': True
#     }]
# peptide = PeptideGenerator.create_peptide(monomers)
# print(peptide)

# template = Chem.MolFromSmiles('CC(C)(OC(OC/C=C/c1cc(C[C@H](C(ON2C(CCC2=O)=O)=O)CC=O)c(F)cc1)=O)C')
# peptide = Chem.MolFromSmiles('N[C@H](Cc1c[nH]c2ccccc12)C(=O)N[C@H](Cc1c[nH]c2ccccc12)C(=O)N[C@H](Cc1ccc(O)cc1)C(=O)O')
# generator = TPHybridGenerator()
# tp_hybrid = generator.connect_template_peptide(generator.modify_template(template), peptide, 'alpha')
# print(tp_hybrid)
