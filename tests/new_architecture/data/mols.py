from rdkit import Chem
from copy import deepcopy

TEST_BACKBONE_1 = {'binary': Chem.MolFromSmiles('N[CH2:1]C(=O)O').ToBinary(
), 'kekule': 'NCC(=O)O', 'mapped_kekule': 'N[CH2:1]C(=O)O'}
TEST_BACKBONE_2 = {'binary': Chem.MolFromSmiles('N[CH2:1]CC(=O)O').ToBinary(
), 'kekule': 'NCCC(=O)O', 'mapped_kekule': 'N[CH2:1]CC(=O)O'}
TEST_BACKBONE_3 = {'binary': Chem.MolFromSmiles('NC[CH2:1]C(=O)O').ToBinary(
), 'kekule': 'NCCC(=O)O', 'mapped_kekule': 'NC[CH2:1]C(=O)O'}

TEST_CONNECTION_1 = {'binary': Chem.MolFromSmiles('C').ToBinary(), 'kekule': 'C'}
TEST_CONNECTION_2 = {'binary': Chem.MolFromSmiles('CC').ToBinary(), 'kekule': 'CC'}

TEST_TEMPLATE_1 = {'binary': Chem.MolFromSmiles(
    'C/C=C/C1=CC(CC[CH:1]=O)=CC=C1').ToBinary(), 'kekule': 'CC(C)(C)OC(=O)OC/C=C/C1=CC(CCC(=O)ON2C(=O)CCC2=O)=CC=C1', 'mapped_kekule': 'C/C=C/C1=CC(CC[CH:1]=O)=CC=C1'}
TEST_TEMPLATE_2 = {'binary': Chem.MolFromSmiles(
    'C/C=C/C1=CC=C(F)C(C[C@@H](CC=O)[CH:1]=O)=C1').ToBinary(), 'kekule': 'CC(C)(C)OC(=O)OC/C=C/C1=CC(C[C@@H](CC=O)C(=O)ON2C(=O)CCC2=O)=C(F)C=C1', 'mapped_kekule': 'C/C=C/C1=CC=C(F)C(C[C@@H](CC=O)[CH:1]=O)=C1'}
TEST_TEMPLATE_3 = {'binary': Chem.MolFromSmiles(
    'C#CCCC[C@](C=O)(CC1=CC=CC(/C=C/C)=C1)C[CH:1]=O').ToBinary(), 'kekule': 'C#CCCC[C@@](Cc1cc(/C=C/COC(OC(C)(C)C)=O)ccc1)(C=O)CC(ON2C(CCC2=O)=O)=O', 'mapped_kekule': 'C#CCCC[C@](C=O)(CC1=CC=CC(/C=C/C)=C1)C[CH:1]=O'}


TEST_SIDECHAIN_1 = {'binary': Chem.MolFromSmiles('CC1=CC=C(O)C=C1').ToBinary(
), 'kekule': 'CC1=CC=C(O)C=C1', 'connection': 'methyl', 'shared_id': 'a'}
TEST_SIDECHAIN_2 = {'binary': Chem.MolFromSmiles('CC1=COC2=NC(=O)[NH]C=C12').ToBinary(
), 'kekule': 'CC1=COC2=NC(=O)[NH]C=C12', 'connection': 'methyl', 'shared_id': 's'}
TEST_SIDECHAIN_3 = {'binary': Chem.MolFromSmiles('CC1=CC(=O)C2=C([NH]1)SC=C2').ToBinary(
), 'kekule': 'CC1=CC(=O)C2=C([NH]1)SC=C2', 'connection': 'methyl', 'shared_id': 'q'}
TEST_SIDECHAIN_4 = {'binary': Chem.MolFromSmiles('CC1=CC=CC2=N[NH]C(=O)N12').ToBinary(
), 'kekule': 'CC1=CC=CC2=N[NH]C(=O)N12', 'connection': 'methyl', 'shared_id': 'af'}

TEST_MONOMER_1 = {'binary': Chem.MolFromSmiles('O=C(O)[C@@H]1C[C@H](OC2=CC=CC=C2)CN1').ToBinary(
), 'kekule': 'O=C(O)[C@@H]1C[C@H](OC2=CC=CC=C2)CN1', 'index': 17, 'backbone': {'_id': 'alpha', 'kekule': 'NCC(=O)O'},
    'sidechain': None, 'connection': None, 'imported': True}
TEST_MONOMER_2 = {'binary': Chem.MolFromSmiles('COC1=CC=C2C(O[C@@H]3CN[C@H](C(=O)O)C3)=CC(C3=CC=CC=C3)=NC2=C1').ToBinary(
), 'kekule': 'COC1=CC=C2C(O[C@@H]3CN[C@H](C(=O)O)C3)=CC(C3=CC=CC=C3)=NC2=C1', 'index': 20, 'backbone': {'_id': 'alpha', 'kekule': 'NCC(=O)O'},
    'sidechain': None, 'connection': None, 'imported': True}
TEST_MONOMER_3 = {'binary': Chem.MolFromSmiles('NC(CC1=CC=CC2=N[NH]C(=O)N12)C(=O)O').ToBinary(
), 'kekule': 'NC(CC1=CC=CC2=N[NH]C(=O)N12)C(=O)O', 'index': 139, 'backbone': {'_id': 'alpha', 'kekule': 'NCC(=O)O'},
    'sidechain': 'af', 'connection': 'methyl', 'imported': False}
TEST_MONOMER_4 = {'binary': Chem.MolFromSmiles('NC(CC(=O)O)CC1=CC=CC2=N[NH]C(=O)N12').ToBinary(
), 'kekule': 'NC(CC(=O)O)CC1=CC=CC2=N[NH]C(=O)N12', 'index': 140, 'backbone': {'_id': 'beta2', 'kekule': 'NCCC(=O)O'},
    'sidechain': 'af', 'connection': 'methyl', 'imported': False}
TEST_MONOMER_5 = {'binary': Chem.MolFromSmiles('NCC(CC1=CC=CC2=N[NH]C(=O)N12)C(=O)O').ToBinary(
), 'kekule': 'NCC(CC1=CC=CC2=N[NH]C(=O)N12)C(=O)O', 'index': 141, 'backbone': {'_id': 'beta3', 'kekule': 'NCCC(=O)O'},
    'sidechain': 'af', 'connection': 'methyl', 'imported': False}
TEST_MONOMER_6 = {'binary': Chem.MolFromSmiles('NCCCC[C@H](N)C(=O)O').ToBinary(), 'kekule': 'NCCCC[C@H](N)C(=O)O',
                  'index': 22, 'backbone': {'_id': 'alpha', 'kekule': 'NCC(=O)O'}, 'sidechain': None,
                  'connection': None, 'imported': True}

TEST_PEPTIDE_1 = {'binary': Chem.MolFromSmiles(
    'NCC(=O)NC(CC(=O)NC(CC(=O)NCC(CC1=CC2=C(OC=C2)S1)C(=O)NC(CC1=CC=C2C=CC=CC=C21)C(=O)O)CC1=CC=CO1)CC1=C2C=CSC2=NS1').ToBinary(),
    'kekule': 'NCC(=O)NC(CC(=O)NC(CC(=O)NCC(CC1=CC2=C(OC=C2)S1)C(=O)NC(CC1=CC=C2C=CC=CC=C21)C(=O)O)CC1=CC=CO1)CC1=C2C=CSC2=NS1',
    'has_cap': False,
    'monomers': [
        {'_id': '12898afefgfad', 'sidechain': 'adwi8', 'proline': False},
        {'_id': 'awfseg4', 'sidechain': '3gdfbv', 'proline': False},
        {'_id': 'asfg43', 'sidechain': 'dws2', 'proline': True}
]}
TEST_PEPTIDE_2 = {'binary': Chem.MolFromSmiles('NC(CC1=CC=CC2=N[NH]C(=O)N12)C(=O)NC(CC(=O)NCC(CC1=CC=CC2=N[NH]C(=O)N12)C(=O)O)CC1=CC=CC2=N[NH]C(=O)N12').ToBinary(
), 'kekule': 'NC(CC1=CC=CC2=N[NH]C(=O)N12)C(=O)NC(CC(=O)NCC(CC1=CC=CC2=N[NH]C(=O)N12)C(=O)O)CC1=CC=CC2=N[NH]C(=O)N12',
    'has_cap': False,
    'monomers': [
    {'_id': '98asfh', 'sidechain': 'af', 'proline': False},
    {'_id': 'adjha82', 'sidechain': 'af', 'proline': False},
    {'_id': 'admaiof7', 'sidechain': 'af', 'proline': False}
]}
TEST_PEPTIDE_3 = {'binary': Chem.MolFromSmiles('COC1=CC=C2C(O[C@H]3C[C@@H](C(=O)O)N(C(=O)C(CC4=CC=CC5=N[NH]C(=O)N45)NC(=O)[C@@H]4C[C@H](OC5=CC(C6=CC=CC=C6)=NC6=CC(OC)=CC=C56)CN4C(=O)[C@@H]4C[C@H](OC5=CC=CC=C5)CN4)C3)=CC(C3=CC=CC=C3)=NC2=C1').ToBinary(
), 'kekule': 'COC1=CC=C2C(O[C@H]3C[C@@H](C(=O)O)N(C(=O)C(CC4=CC=CC5=N[NH]C(=O)N45)NC(=O)[C@@H]4C[C@H](OC5=CC(C6=CC=CC=C6)=NC6=CC(OC)=CC=C56)CN4C(=O)[C@@H]4C[C@H](OC5=CC=CC=C5)CN4)C3)=CC(C3=CC=CC=C3)=NC2=C1',
    'has_cap': False,
    'monomers': [
    {'_id': 'ad98fh', 'sidechain': None, 'proline': True},
    {'_id': 'sdwd89cvh', 'sidechain': None, 'proline': True},
    {'_id': '98asfh', 'sidechain': 'af', 'proline': False},
    {'_id': 'sdwd89cvh', 'sidechain': None, 'proline': True}
]}
TEST_PEPTIDE_4_CAP = {'binary': Chem.MolFromSmiles('COC1=CC=C2C(O[C@H]3C[C@@H](C(=O)NCCC4=CC=CC5=N[NH]C(=O)N45)N(C(=O)C(CC4=CC=CC5=N[NH]C(=O)N45)NC(=O)[C@@H]4C[C@H](OC5=CC(C6=CC=CC=C6)=NC6=CC(OC)=CC=C56)CN4C(=O)[C@@H]4C[C@H](OC5=CC=CC=C5)CN4)C3)=CC(C3=CC=CC=C3)=NC2=C1').ToBinary(
), 'kekule': 'COC1=CC=C2C(O[C@H]3C[C@@H](C(=O)NCCC4=CC=CC5=N[NH]C(=O)N45)N(C(=O)C(CC4=CC=CC5=N[NH]C(=O)N45)NC(=O)[C@@H]4C[C@H](OC5=CC(C6=CC=CC=C6)=NC6=CC(OC)=CC=C56)CN4C(=O)[C@@H]4C[C@H](OC5=CC=CC=C5)CN4)C3)=CC(C3=CC=CC=C3)=NC2=C1',
    'has_cap': True,
    'monomers': [
    {'_id': 'ad98fh', 'sidechain': None, 'proline': True},
    {'_id': 'sdwd89cvh', 'sidechain': None, 'proline': True},
    {'_id': '98asfh', 'sidechain': 'af', 'proline': False},
    {'_id': 'sdwd89cvh', 'sidechain': None, 'proline': True},
    {'_id': '98asfh', 'sidechain': 'af', 'proline': False}
]}
TEST_PEPTIDE_4 = {'binary': Chem.MolFromSmiles('COC1=CC=C2C(O[C@H]3C[C@@H](C(=O)NC(CC4=CC=CC5=N[NH]C(=O)N45)C(=O)O)N(C(=O)C(CC4=CC=CC5=N[NH]C(=O)N45)NC(=O)[C@@H]4C[C@H](OC5=CC(C6=CC=CC=C6)=NC6=CC(OC)=CC=C56)CN4C(=O)[C@@H]4C[C@H](OC5=CC=CC=C5)CN4)C3)=CC(C3=CC=CC=C3)=NC2=C1').ToBinary(
), 'kekule': 'COC1=CC=C2C(O[C@H]3C[C@@H](C(=O)NC(CC4=CC=CC5=N[NH]C(=O)N45)C(=O)O)N(C(=O)C(CC4=CC=CC5=N[NH]C(=O)N45)NC(=O)[C@@H]4C[C@H](OC5=CC(C6=CC=CC=C6)=NC6=CC(OC)=CC=C56)CN4C(=O)[C@@H]4C[C@H](OC5=CC=CC=C5)CN4)C3)=CC(C3=CC=CC=C3)=NC2=C1',
    'has_cap': False,
    'monomers': [
    {'_id': 'ad98fh', 'sidechain': None, 'proline': True},
    {'_id': 'sdwd89cvh', 'sidechain': None, 'proline': True},
    {'_id': '98asfh', 'sidechain': 'af', 'proline': False},
    {'_id': 'sdwd89cvh', 'sidechain': None, 'proline': True},
    {'_id': '98asfh', 'sidechain': 'af', 'proline': False}
]}
TEST_PEPTIDE_5 = {'binary': Chem.MolFromSmiles('COC1=CC=C2C(O[C@H]3C[C@@H](C(=O)N[C@@H](CCCCN)C(=O)O)N(C(=O)C(CC4=CC=CC5=N[NH]C(=O)N45)NC(=O)[C@@H]4C[C@H](OC5=CC(C6=CC=CC=C6)=NC6=CC(OC)=CC=C56)CN4C(=O)[C@@H]4C[C@H](OC5=CC=CC=C5)CN4)C3)=CC(C3=CC=CC=C3)=NC2=C1').ToBinary(
), 'kekule': 'COC1=CC=C2C(O[C@H]3C[C@@H](C(=O)N[C@@H](CCCCN)C(=O)O)N(C(=O)C(CC4=CC=CC5=N[NH]C(=O)N45)NC(=O)[C@@H]4C[C@H](OC5=CC(C6=CC=CC=C6)=NC6=CC(OC)=CC=C56)CN4C(=O)[C@@H]4C[C@H](OC5=CC=CC=C5)CN4)C3)=CC(C3=CC=CC=C3)=NC2=C1',
    'has_cap': False,
    'monomers': [
    {'_id': 'ad98fh', 'sidechain': None, 'proline': True},
    {'_id': 'sdwd89cvh', 'sidechain': None, 'proline': True},
    {'_id': '98asfh', 'sidechain': 'af', 'proline': False},
    {'_id': 'sdwd89cvh', 'sidechain': None, 'proline': True},
    {'_id': '12fakwd', 'sidechain': None, 'proline': False}
]}


TEST_PEPTIDE_1_WITH_ID = deepcopy(TEST_PEPTIDE_1)
TEST_PEPTIDE_5_WITH_ID = deepcopy(TEST_PEPTIDE_5)
TEST_PEPTIDE_1_WITH_ID['_id'] = 'aefoi249'
TEST_PEPTIDE_5_WITH_ID['_id'] = 'fasef00'
TEST_PEPTIDE_1_WITH_ID.pop('binary')
TEST_PEPTIDE_5_WITH_ID.pop('binary')
TEST_TEMPLATE_PEPTIDE_1 = {'binary': Chem.MolFromSmiles('C/C=C/C1=CC=CC(CCC(=O)NC(CCC2=C3C=COC3=CS2)C(=O)NC(CCC2=CC=CC3=COC=C23)C(=O)NC(CCC2=CC=C(O)C=C2)C(=O)N2[C@H](C(=O)NC(C(=O)O)C(C)C)C[C@H]3C[C@H]32)=C1').ToBinary(),
                           'kekule': 'C/C=C/C1=CC=CC(CCC(=O)NC(CCC2=C3C=COC3=CS2)C(=O)NC(CCC2=CC=CC3=COC=C23)C(=O)NC(CCC2=CC=C(O)C=C2)C(=O)N2[C@H](C(=O)NC(C(=O)O)C(C)C)C[C@H]3C[C@H]32)=C1',
                           'template': 'temp1',
                           'peptide': TEST_PEPTIDE_1_WITH_ID}
TEST_TEMPLATE_PEPTIDE_2 = {'binary': Chem.MolFromSmiles('C/C=C/C1=CC=CC(CCC(=O)N2C[C@@H](OC3=CC=CC=C3)C[C@H]2C(=O)N2C[C@@H](OC3=CC(C4=CC=CC=C4)=NC4=CC(OC)=CC=C34)C[C@H]2C(=O)NC(CC2=CC=CC3=N[NH]C(=O)N23)C(=O)N2C[C@@H](OC3=CC(C4=CC=CC=C4)=NC4=CC(OC)=CC=C34)C[C@H]2C(=O)N[C@@H](CCCCN)C(=O)O)=C1').ToBinary(),
                           'kekule': 'C/C=C/C1=CC=CC(CCC(=O)N2C[C@@H](OC3=CC=CC=C3)C[C@H]2C(=O)N2C[C@@H](OC3=CC(C4=CC=CC=C4)=NC4=CC(OC)=CC=C34)C[C@H]2C(=O)NC(CC2=CC=CC3=N[NH]C(=O)N23)C(=O)N2C[C@@H](OC3=CC(C4=CC=CC=C4)=NC4=CC(OC)=CC=C34)C[C@H]2C(=O)N[C@@H](CCCCN)C(=O)O)=C1',
                           'template': 'temp1',
                           'peptide': TEST_PEPTIDE_5_WITH_ID}
TEST_TEMPLATE_PEPTIDE_3 = {'binary': Chem.MolFromSmiles('C/C=C/C1=CC=CC(CCC(=O)NCCCC[C@H](NC(=O)[C@@H]2C[C@H](OC3=CC(C4=CC=CC=C4)=NC4=CC(OC)=CC=C34)CN2C(=O)C(CC2=CC=CC3=N[NH]C(=O)N23)NC(=O)[C@@H]2C[C@H](OC3=CC(C4=CC=CC=C4)=NC4=CC(OC)=CC=C34)CN2C(=O)[C@@H]2C[C@H](OC3=CC=CC=C3)CN2)C(=O)O)=C1').ToBinary(),
                           'kekule': 'C/C=C/C1=CC=CC(CCC(=O)NCCCC[C@H](NC(=O)[C@@H]2C[C@H](OC3=CC(C4=CC=CC=C4)=NC4=CC(OC)=CC=C34)CN2C(=O)C(CC2=CC=CC3=N[NH]C(=O)N23)NC(=O)[C@@H]2C[C@H](OC3=CC(C4=CC=CC=C4)=NC4=CC(OC)=CC=C34)CN2C(=O)[C@@H]2C[C@H](OC3=CC=CC=C3)CN2)C(=O)O)=C1',
                           'template': 'temp1',
                           'peptide': TEST_PEPTIDE_5_WITH_ID}
TEST_TEMPLATE_PEPTIDE_4 = {'binary': Chem.MolFromSmiles('C/C=C/C1=CC=C(F)C(C[C@@H](CC=O)C(=O)N2C[C@@H](OC3=CC=CC=C3)C[C@H]2C(=O)N2C[C@@H](OC3=CC(C4=CC=CC=C4)=NC4=CC(OC)=CC=C34)C[C@H]2C(=O)NC(CC2=CC=CC3=N[NH]C(=O)N23)C(=O)N2C[C@@H](OC3=CC(C4=CC=CC=C4)=NC4=CC(OC)=CC=C34)C[C@H]2C(=O)N[C@@H](CCCCN)C(=O)O)=C1').ToBinary(),
                           'kekule': 'C/C=C/C1=CC=C(F)C(C[C@@H](CC=O)C(=O)N2C[C@@H](OC3=CC=CC=C3)C[C@H]2C(=O)N2C[C@@H](OC3=CC(C4=CC=CC=C4)=NC4=CC(OC)=CC=C34)C[C@H]2C(=O)NC(CC2=CC=CC3=N[NH]C(=O)N23)C(=O)N2C[C@@H](OC3=CC(C4=CC=CC=C4)=NC4=CC(OC)=CC=C34)C[C@H]2C(=O)N[C@@H](CCCCN)C(=O)O)=C1',
                           'template': 'temp2',
                           'peptide': TEST_PEPTIDE_5_WITH_ID}
TEST_TEMPLATE_PEPTIDE_5 = {'binary': Chem.MolFromSmiles('C/C=C/C1=CC=C(F)C(C[C@@H](CC=O)C(=O)NCCCC[C@H](NC(=O)[C@@H]2C[C@H](OC3=CC(C4=CC=CC=C4)=NC4=CC(OC)=CC=C34)CN2C(=O)C(CC2=CC=CC3=N[NH]C(=O)N23)NC(=O)[C@@H]2C[C@H](OC3=CC(C4=CC=CC=C4)=NC4=CC(OC)=CC=C34)CN2C(=O)[C@@H]2C[C@H](OC3=CC=CC=C3)CN2)C(=O)O)=C1').ToBinary(),
                           'kekule': 'C/C=C/C1=CC=C(F)C(C[C@@H](CC=O)C(=O)NCCCC[C@H](NC(=O)[C@@H]2C[C@H](OC3=CC(C4=CC=CC=C4)=NC4=CC(OC)=CC=C34)CN2C(=O)C(CC2=CC=CC3=N[NH]C(=O)N23)NC(=O)[C@@H]2C[C@H](OC3=CC(C4=CC=CC=C4)=NC4=CC(OC)=CC=C34)CN2C(=O)[C@@H]2C[C@H](OC3=CC=CC=C3)CN2)C(=O)O)=C1',
                           'template': 'temp2',
                           'peptide': TEST_PEPTIDE_5_WITH_ID}
TEST_TEMPLATE_PEPTIDE_6 = {'binary': Chem.MolFromSmiles('C#CCCC[C@@](C=O)(CC(=O)N1C[C@@H](OC2=CC=CC=C2)C[C@H]1C(=O)N1C[C@@H](OC2=CC(C3=CC=CC=C3)=NC3=CC(OC)=CC=C23)C[C@H]1C(=O)NC(CC1=CC=CC2=N[NH]C(=O)N12)C(=O)N1C[C@@H](OC2=CC(C3=CC=CC=C3)=NC3=CC(OC)=CC=C23)C[C@H]1C(=O)N[C@@H](CCCCN)C(=O)O)CC1=CC=CC(/C=C/C)=C1').ToBinary(),
                           'kekule': 'C#CCCC[C@@](C=O)(CC(=O)N1C[C@@H](OC2=CC=CC=C2)C[C@H]1C(=O)N1C[C@@H](OC2=CC(C3=CC=CC=C3)=NC3=CC(OC)=CC=C23)C[C@H]1C(=O)NC(CC1=CC=CC2=N[NH]C(=O)N12)C(=O)N1C[C@@H](OC2=CC(C3=CC=CC=C3)=NC3=CC(OC)=CC=C23)C[C@H]1C(=O)N[C@@H](CCCCN)C(=O)O)CC1=CC=CC(/C=C/C)=C1',
                           'template': 'temp3',
                           'peptide': TEST_PEPTIDE_5_WITH_ID}
TEST_TEMPLATE_PEPTIDE_7 = {'binary': Chem.MolFromSmiles('C#CCCC[C@@](C=O)(CC(=O)NCCCC[C@H](NC(=O)[C@@H]1C[C@H](OC2=CC(C3=CC=CC=C3)=NC3=CC(OC)=CC=C23)CN1C(=O)C(CC1=CC=CC2=N[NH]C(=O)N12)NC(=O)[C@@H]1C[C@H](OC2=CC(C3=CC=CC=C3)=NC3=CC(OC)=CC=C23)CN1C(=O)[C@@H]1C[C@H](OC2=CC=CC=C2)CN1)C(=O)O)CC1=CC=CC(/C=C/C)=C1').ToBinary(),
                           'kekule': 'C#CCCC[C@@](C=O)(CC(=O)NCCCC[C@H](NC(=O)[C@@H]1C[C@H](OC2=CC(C3=CC=CC=C3)=NC3=CC(OC)=CC=C23)CN1C(=O)C(CC1=CC=CC2=N[NH]C(=O)N12)NC(=O)[C@@H]1C[C@H](OC2=CC(C3=CC=CC=C3)=NC3=CC(OC)=CC=C23)CN1C(=O)[C@@H]1C[C@H](OC2=CC=CC=C2)CN1)C(=O)O)CC1=CC=CC(/C=C/C)=C1',
                           'template': 'temp3',
                           'peptide': TEST_PEPTIDE_5_WITH_ID}
