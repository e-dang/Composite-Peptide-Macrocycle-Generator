from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from utils import Database, generate_rxn_temp, generate_candidates, write_mols, read_mols

mol = Chem.MolFromSmiles('CC(C)(C)OC(=O)OC/C=C/c1ccc(F)c(C[C@@H](CC2OCCCN2C(=O)OC(C)(C)C)C(=O)O)c1')
# mol = Chem.MolFromSmarts('N[C@H](Cc1cc2ccc(=O)[nH]c2s1)C(=O)N1CCC[C@H:2]1C=O')
# Draw.MolToImage(mol, size=(1000, 1000), includeAtomNumbers=True).show()
Draw.MolToImage(mol, size=(500, 500)).show()


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
