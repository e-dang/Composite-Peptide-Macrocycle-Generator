from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from utils import Database, write_mols, read_mols
import json


# rxn = AllChem.ReactionFromSmarts(
# '([*:3]Cc1cccc(/C=C/[C:1]OC(=O)OC(C)(C)C)c1.c1oc(C[*:4])c[cH:2]1)>>C(=C/[C:1][c:2]1coc(C[*:4])c1)\\c1cccc(C[*:3])c1')
# mol = Chem.MolFromSmiles('[CH3][CH2]')
# mol = Chem.MolFromSmiles('O=c1nc2c(co1)c(C[*:4])c[nH:2]2')
# Draw.MolToImage(mol, size=(1000, 1000), includeAtomNumbers=True).show()
# Draw.MolToImage(mol, size=(500, 500)).show()
# Draw.ReactionToImage(rxn, subImgSize=(500, 500)).show()


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

# with open('/Users/ericdang/Documents/UCLA_Research/macrocycles/smiles/pre_monomers/', 'w') as fout:
#     mols = read_mols('/Users/ericdang/Documents/UCLA_Research/macrocycles/chemdraw/monomer_backbone.sdf')
#     collection = []
#     for mol in mols:
#         doc = {}
#         doc['type'] = 'alpha_amino_acid'
#         doc['smiles'] = Chem.MolToSmiles(mol)
#         collection.append(doc)
#     json.dump(collection, fout)


# mols = read_mols('/Users/ericdang/Documents/UCLA_Research/macrocycles/chemdraw/D-amino_acids.sdf')
# for mol in mols:
#     print(Chem.MolToSmiles(mol))

mols = read_mols('/Users/ericdang/Documents/UCLA_Research/macrocycles/chemdraw/amino_acid_side_chains.sdf')
for mol in mols:
    print(Chem.MolToSmiles(mol))
    Draw.MolToImage(mol).show()
    l = input()
