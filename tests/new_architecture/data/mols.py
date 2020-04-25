from rdkit import Chem
from copy import deepcopy

TEST_BACKBONE_1 = {'binary': Chem.MolFromSmiles('N[CH2:1]C(=O)O').ToBinary(
), 'kekule': 'N[CH2:1]C(=O)O'}
TEST_BACKBONE_2 = {'binary': Chem.MolFromSmiles('N[CH2:1]CC(=O)O').ToBinary(
), 'kekule': 'N[CH2:1]CC(=O)O'}
TEST_BACKBONE_3 = {'binary': Chem.MolFromSmiles('NC[CH2:1]C(=O)O').ToBinary(
), 'kekule': 'NC[CH2:1]C(=O)O'}

TEST_CONNECTION_1 = {'binary': Chem.MolFromSmiles('CC').ToBinary(), 'kekule': 'CC'}

TEST_TEMPLATE_1 = {'binary': Chem.MolFromSmiles(
    'CC(C)(C)OC(=O)OC/C=C/C1=CC(CCC(=O)ON2C(=O)CCC2=O)=CC=C1').ToBinary(), 'kekule': 'CC(C)(C)OC(=O)OC/C=C/C1=CC(CCC(=O)ON2C(=O)CCC2=O)=CC=C1'}
TEST_TEMPLATE_2 = {'binary': Chem.MolFromSmiles(
    'CC(C)(C)OC(=O)OC/C=C/C1=CC(C[C@@H](CC=O)C(=O)ON2C(=O)CCC2=O)=C(F)C=C1').ToBinary(), 'kekule': 'CC(C)(C)OC(=O)OC/C=C/C1=CC(C[C@@H](CC=O)C(=O)ON2C(=O)CCC2=O)=C(F)C=C1'}
TEST_TEMPLATE_3 = {'binary': Chem.MolFromSmiles(
    'C#CCCC[C@@](Cc1cc(/C=C/COC(OC(C)(C)C)=O)ccc1)(C=O)CC(ON2C(CCC2=O)=O)=O').ToBinary(), 'kekule': 'C#CCCC[C@@](Cc1cc(/C=C/COC(OC(C)(C)C)=O)ccc1)(C=O)CC(ON2C(CCC2=O)=O)=O'}


TEST_SIDECHAIN_1 = {'binary': Chem.MolFromSmiles('CC1=CC=C(O)C=C1').ToBinary(
), 'kekule': 'CC1=CC=C(O)C=C1', 'connection': 'methyl', 'shared_id': 'a'}
TEST_SIDECHAIN_2 = {'binary': Chem.MolFromSmiles('CC1=COC2=NC(=O)[NH]C=C12').ToBinary(
), 'kekule': 'CC1=COC2=NC(=O)[NH]C=C12', 'connection': 'methyl', 'shared_id': 's'}
TEST_SIDECHAIN_3 = {'binary': Chem.MolFromSmiles('CC1=CC(=O)C2=C([NH]1)SC=C2').ToBinary(
), 'kekule': 'CC1=CC(=O)C2=C([NH]1)SC=C2', 'connection': 'methyl', 'shared_id': 'q'}


TEST_MONOMER_1 = {'binary': Chem.MolFromSmiles('O=C(O)[C@@H]1C[C@H](OC2=CC=CC=C2)CN1').ToBinary(
), 'kekule': 'O=C(O)[C@@H]1C[C@H](OC2=CC=CC=C2)CN1', 'index': 17, 'required': True, 'backbone': 'alpha',
    'sidechain': None, 'connection': None, 'is_proline': True, 'imported': True}
TEST_MONOMER_2 = {'binary': Chem.MolFromSmiles('COC1=CC=C2C(O[C@@H]3CN[C@H](C(=O)O)C3)=CC(C3=CC=CC=C3)=NC2=C1').ToBinary(
), 'kekule': 'COC1=CC=C2C(O[C@@H]3CN[C@H](C(=O)O)C3)=CC(C3=CC=CC=C3)=NC2=C1', 'index': 20, 'required': True, 'backbone': 'alpha',
    'sidechain': None, 'connection': None, 'is_proline': True, 'imported': True}
TEST_MONOMER_3 = {'binary': Chem.MolFromSmiles('NCC(CC1=CC=CC2=N[NH]C(=O)N12)C(=O)O').ToBinary(
), 'kekule': 'NCC(CC1=CC=CC2=N[NH]C(=O)N12)C(=O)O', 'index': 140, 'required': True, 'backbone': 'beta3',
    'sidechain': 'af', 'connection': 'methyl', 'is_proline': False, 'imported': False}


TEST_PEPTIDE_1 = {'binary': Chem.MolFromSmiles(
    'NCC(=O)NC(CC(=O)NC(CC(=O)NCC(CC1=CC2=C(OC=C2)S1)C(=O)NC(CC1=CC=C2C=CC=CC=C21)C(=O)O)CC1=CC=CO1)CC1=C2C=CSC2=NS1').ToBinary(),
    'kekule': 'NCC(=O)NC(CC(=O)NC(CC(=O)NCC(CC1=CC2=C(OC=C2)S1)C(=O)NC(CC1=CC=C2C=CC=CC=C21)C(=O)O)CC1=CC=CO1)CC1=C2C=CSC2=NS1',
    'has_cap': False,
    'monomers': [
        {'_id': '12898afefgfad', 'sidechain': 'adwi8', 'is_proline': False},
        {'_id': 'awfseg4', 'sidechain': '3gdfbv', 'is_proline': False},
        {'_id': 'asfg43', 'sidechain': 'dws2', 'is_proline': True}
]}

TEST_PEPTIDE_WITH_ID = deepcopy(TEST_PEPTIDE_1)
TEST_PEPTIDE_WITH_ID['_id'] = 'aefoi249'
TEST_PEPTIDE_WITH_ID.pop('binary')
TEST_TEMPLATE_PEPTIDE_1 = {'binary': Chem.MolFromSmiles('C/C=C/C1=CC=CC(CCC(=O)NC(CCC2=C3C=COC3=CS2)C(=O)NC(CCC2=CC=CC3=COC=C23)C(=O)NC(CCC2=CC=C(O)C=C2)C(=O)N2[C@H](C(=O)NC(C(=O)O)C(C)C)C[C@H]3C[C@H]32)=C1').ToBinary(),
                           'kekule': 'C/C=C/C1=CC=CC(CCC(=O)NC(CCC2=C3C=COC3=CS2)C(=O)NC(CCC2=CC=CC3=COC=C23)C(=O)NC(CCC2=CC=C(O)C=C2)C(=O)N2[C@H](C(=O)NC(C(=O)O)C(C)C)C[C@H]3C[C@H]32)=C1',
                           'template': 'temp1',
                           'peptide': TEST_PEPTIDE_WITH_ID}
