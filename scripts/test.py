import sys
# from rdkit import Chem
# from rdkit.Chem import AllChem
# from pprint import pprint
# from rdkit.Chem import Draw
# from rdkit.Chem import AllChem
# # from utils.database import collection_example_reactions_details, collection_templates, collection_candidates, collection_example_reactions_smilesonly
# import os
# from indigo.indigo import *
# from indigo.indigo_renderer import *
# import re
# from utils import Database, generate_rxn_temp, generate_candidates, write_mols, read_mols

# indigo = Indigo()


# rxn = indigo.createReaction()
# reactant = indigo.loadMolecule('CC(C)[C@H](NC(=O)CCc1cccc(/C=C/COC(=O)OC(C)(C)C)c1)C(=O)N[C@@H](CCC(N)=O)C(=O)N[C@@H](Cc1ccc(O)cc1)C(=O)N[C@@H](Cc1c[nH]cn1)C(N)=O')
# product = indigo.loadMolecule('CC(C)[C@@H]1NC(=O)CCc2cccc(c2)/C=C/Cn2cnc(c2)C[C@@H](C(N)=O)NC(=O)[C@H](Cc2ccc(O)cc2)NC(=O)[C@H](CCC(N)=O)NC1=O')
# rxn.addReactant(reactant)
# rxn.addProduct(product)
# rxn.automap('discard')
#
# # print(Indigo.RC_NOT_CENTER) #-1
# # print(Indigo.RC_UNMARKED)   #0
# # print(Indigo.RC_CENTER) #1
# # print(Indigo.RC_UNCHANGED)  #2
# # print(Indigo.RC_MADE_OR_BROKEN) #4
# # print(Indigo.RC_ORDER_CHANGED)  #8
#
# rxn.correctReactingCenters()
# round = 0
# for m in rxn.iterateMolecules():
#     if round == 0:
#         round += 1
#         continue
#
#     for b in m.iterateBonds():
#         print(rxn.reactingCenter(b))
#         if rxn.reactingCenter(b) == 4:
#             atom = b.destination()
#             # if atom.symbol() == 'C':
#             #     for nei in atom.iterateNeighbors():
#             #         print(nei.symbol(), rxn.atomMappingNumber(nei))
#             # exit()
#             print(b.source().symbol(), rxn.atomMappingNumber(b.source()))
#             print(b.destination().symbol(), rxn.atomMappingNumber(b.destination()))





# mol = Chem.MolFromSmiles('C[C@@H](O)[C@H](NC(=O)[C@H](Cc1scc2[nH]cnc12)NC(=O)[C@H](CCCNC(=N)N)NC(=O)CCc1cccc(/C=C/COC(=O)OC(C)(C)C)c1)C(=O)N[C@@H](Cc1cc2[nH]ccc(=O)c2s1)C(=O)N[C@@H](CCC(N)=O)C(=O)O')
# mol = Chem.MolFromSmarts('c1c[cH:2]oc1C[*:4]')
# print(mol.HasSubstructMatch(patt, useChirality=False))
# Draw.MolToImage(mol, size=(1000,1000), includeAtomNumbers=True).show()
# Draw.MolToImage(mol, size=(500,500)).show()
# Chem.MolFromSmarts('[C,C(*)]')

# data = {'name': 'p-amino-phe', 'smiles': '[*:4]Cc1ccc([N:2])cc1', 'chain_ind': '4', 'rxn_ind': '2'}
# db = Database()
# temp = db.get_doc('templates', {'name': 'temp1a'})
# side_chain = db.get_doc('side_chains', {})
# for sc in side_chain:
#     print(sc['name'])
#     print(generate_rxn_temp(temp[0], sc, store=True))

# db.insert('side_chain', data)

# Draw.ReactionToImage(AllChem.ReactionFromSmarts('([*:3]Cc1cccc(/C=C/[C:1]OC(=O)OC(C)(C)C)c1.[*:4]Cc1ccc([O:2])cc1)>>([*:4]Cc1ccc([O:2][C:1]/C=C/c2cccc(C[*:3])c2)cc1)'), subImgSize=(500,500)).show()
# rxn = AllChem.ReactionFromSmarts('([*:3]Cc1cccc(/C=C/[C:1]OC(=O)OC(C)(C)C)c1.[*:4]Cc1ccc([O:2])cc1)>>([*:4]Cc1ccc([O:2][C:1]/C=C/c2cccc(C[*:3])c2)cc1)')
# prod = rxn.RunReactants((Chem.MolFromSmiles('CC(C)[C@H](NC(=O)CCc1cccc(/C=C/COC(=O)OC(C)(C)C)c1)C(=O)N[C@@H](CCC(N)=O)C(=O)N[C@@H](Cc1ccc(O)cc1)C(=O)N[C@@H](Cc1c[nH]cn1)C(N)=O'),))
# # #
# print(len(prod))
# unique = {}
# for x in prod:
#     # print(type(x[0]))
#     Chem.SanitizeMol(x[0])
#     Draw.MolToImage(x[0]).show()
#     smi = Chem.MolToSmiles(x[0])
#     unique[smi] = x[0]
#
# print(sorted(unique.keys()))

# reactants = read_mols(verbose=True)
# for reactant in reactants:
#     print(Chem.MolToSmiles(reactant))

# reactant = 'CSCC[C@H](NC(=O)[C@H](CCCNC(=O)CCc1cccc(/C=C/COC(=O)OC(C)(C)C)c1)CC(=O)[C@@H](N)C(C)C)C(=O)N[C@@H](Cc1ccc(O)cc1)C(N)=O'
# for reactant in reactants:
#     reactant = Chem.MolToSmiles(reactant)
#     products = generate_candidates(reactant)
#     writer = Chem.SDWriter('/Users/ericdang/Documents/UCLA_Research/chemdraw/candidates.sdf')
#     writer.write(Chem.MolFromSmiles(reactant))
#     for mol in products:
#         writer.write(Chem.MolFromSmiles(mol))
#     print(products)


# db = Database()
# for ind, cand in enumerate(db.find_all('candidates')):
#     if ind >= 12:
#         write_mols(cand['products'], 'candidates/demo_' + str(ind) + '.sdf')


print(sys.version)
print(sys.executable)
print


