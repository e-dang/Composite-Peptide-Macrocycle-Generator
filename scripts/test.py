from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import rdmolfiles
from utils import Database, write_mols, read_mols, create_monomer_requirements
import json


rxn = AllChem.ReactionFromSmarts(
    '([*:3]Cc1cccc(/C=C/[C:1]OC(=O)OC(C)(C)C)c1.c1oc[cH:2]c1C[*:4])>>C(=C/[C:1][c:2]1cocc1C[*:4])\\c1cccc(C[*:3])c1')
# mol = Chem.MolFromSmarts('c1nc2n(n1)c(CC[CH3:4][*:4])n[nH:2]2')
# mol = Chem.MolFromSmiles('Cc1ccc(O)cc1')
# Draw.MolToImage(mol, size=(1000, 1000), includeAtomNumbers=True).show()
# Draw.MolToImage(mol, size=(500, 500)).show()
Draw.ReactionToImage(rxn, subImgSize=(500, 500)).show()

# mol = Chem.MolFromSmiles(
# 'O=C1CCc2cccc(c2)/C=C/Cc2cccc3[nH]cc(c23)C[C@@H](C(=O)N[C@@H](Cc2ccc(O)cc2)C(=O)O)NC(=O)[C@H](Cc2c[nH]c3ccccc23)N1')
# print(Chem.MolToSmiles(mol))
# tmp = AllChem.Compute2DCoords(mol)
# print(tmp)
# AllChem.Ge
# AllChem.AddHs(mol)
# AllChem.EmbedMolecule(mol)
# print(AllChem.MMFFOptimizeMolecule(mol, maxIters=1000))
# rdmolfiles.MolToPDBFile(mol, 'test.pdb')
# writer = Chem.SDWriter('test.sdf')
# writer.write(mol)
# Draw.
# Draw.MolToImage(mol).show()

# reactants = read_mols(verbose=False)
# for reactant in reactants:
#     print(Chem.MolToSmiles(reactant))

# aa = "N=C(N)NCCC[C@H](N)C(=O)O, NC(=O)C[C@H](N)C(=O)O, N[C@@H](CC(=O)O)C(=O)O, N[C@@H](CS)C(=O)O, N[C@@H](CCC(=O)O)C(=O)O, C[C@H](N)C(=O)O, NCC(=O)O, N[C@@H](Cc1c[nH]cn1)C(=O)O, CC[C@H](C)[C@H](N)C(=O)O, CC(C)C[C@H](N)C(=O)O, NCCCC[C@H](N)C(=O)O, CSCC[C@H](N)C(=O)O, N[C@@H](Cc1ccccc1)C(=O)O, O=C(O)[C@@H]1CCCN1, N[C@@H](CO)C(=O)O, C[C@@H](O)[C@H](N)C(=O)O, N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O, N[C@@H](Cc1ccc(O)cc1)C(=O)O, CC(C)[C@H](N)C(=O)O".split(', ')
# for x in aa:
#     Draw.MolToImage(Chem.MolFromSmiles(x)).show()

# db = Database()
# for ind, cand in enumerate(db.find_all('candidates')):
#     if ind >= 12:
#         write_mols(cand['products'], 'candidates/demo_' + str(ind) + '.sdf')

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


# mols = read_mols('/Users/ericdang/Documents/UCLA_Research/macrocycles/chemdraw/D-amino_acids.sdf')
# for mol in mols:
#     print(Chem.MolToSmiles(mol))

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

# mols = read_mols('/Users/ericdang/Documents/UCLA_Research/macrocycles/chemdraw/pre_monomer/side_chains_likely1.sdf')
# d = {}
# for mol in mols:
#     smiles = Chem.MolToSmiles(mol)
#     if smiles not in d.keys():
#         d[smiles] = 1
#     else:
#         d[smiles] += 1

# for key, value in d.items():
#     if value >= 2:
#         print(key, value)

# create_monomer_requirements()
