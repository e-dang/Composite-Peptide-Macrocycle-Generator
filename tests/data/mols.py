import uuid
from copy import deepcopy

from rdkit import Chem
from rdkit.Chem import AllChem

import cpmg.models as models
import cpmg.reactions as rxns

TEST_BACKBONE_1 = {'binary': Chem.MolFromSmiles('N[CH2:1]C(=O)O').ToBinary(
), 'kekule': 'NCC(=O)O', 'mapped_kekule': 'N[CH2:1]C(=O)O'}
TEST_BACKBONE_2 = {'binary': Chem.MolFromSmiles('N[CH2:1]CC(=O)O').ToBinary(
), 'kekule': 'NCCC(=O)O', 'mapped_kekule': 'N[CH2:1]CC(=O)O'}
TEST_BACKBONE_3 = {'binary': Chem.MolFromSmiles('NC[CH2:1]C(=O)O').ToBinary(
), 'kekule': 'NCCC(=O)O', 'mapped_kekule': 'NC[CH2:1]C(=O)O'}
TEST_BACKBONES = [TEST_BACKBONE_1, TEST_BACKBONE_2, TEST_BACKBONE_3]

TEST_CONNECTION_1 = {'binary': Chem.MolFromSmiles('C').ToBinary(), 'kekule': 'C'}
TEST_CONNECTION_2 = {'binary': Chem.MolFromSmiles('CC').ToBinary(), 'kekule': 'CC'}
TEST_CONNECTIONS = [TEST_CONNECTION_1, TEST_CONNECTION_2]

TEST_TEMPLATE_1 = {'binary': Chem.MolFromSmiles(
    'CC(C)(C)OC(=O)OC/C=C/C1=CC(CCC(=O)ON2C(=O)CCC2=O)=CC=C1').ToBinary(),
    'kekule': 'CC(C)(C)OC(=O)OC/C=C/C1=CC(CCC(=O)ON2C(=O)CCC2=O)=CC=C1',
    'oligomerization_kekule': f'C/C=C/C1=CC(CC[CH:{models.Template.OLIGOMERIZATION_MAP_NUM}]=O)=CC=C1',
    'friedel_crafts_kekule': f'[*:{models.Template.WC_MAP_NUM_1}]/C=C/[CH3:{models.Template.EAS_MAP_NUM}]',
    'tsuji_trost_kekule': f'[*:{models.Template.WC_MAP_NUM_1}]/C=C/[CH3:{models.Template.EAS_MAP_NUM}]',
    'pictet_spangler_kekule': None,
    'template_pictet_spangler_kekule': None,
    'pyrroloindoline_kekule': f'[*:{models.Template.WC_MAP_NUM_1}]/C=C/[CH3:{models.Template.EAS_MAP_NUM}]',
    'aldehyde_cyclization_kekule': None}
TEST_TEMPLATE_2 = {'binary': Chem.MolFromSmiles(
    'CC(C)(C)OC(=O)OC/C=C/C1=CC(C[C@@H](CC=O)C(=O)ON2C(=O)CCC2=O)=C(F)C=C1').ToBinary(),
    'kekule': 'CC(C)(C)OC(=O)OC/C=C/C1=CC(C[C@@H](CC=O)C(=O)ON2C(=O)CCC2=O)=C(F)C=C1',
    'oligomerization_kekule': f'C/C=C/C1=CC=C(F)C(C[C@@H](CC=O)[CH:{models.Template.OLIGOMERIZATION_MAP_NUM}]=O)=C1',
    'friedel_crafts_kekule': f'[*:{models.Template.WC_MAP_NUM_1}]/C=C/[CH3:{models.Template.EAS_MAP_NUM}]',
    'tsuji_trost_kekule': f'[*:{models.Template.WC_MAP_NUM_1}]/C=C/[CH3:{models.Template.EAS_MAP_NUM}]',
    'pictet_spangler_kekule': f'[O:{models.Template.PS_OXYGEN_MAP_NUM}]=[CH1:{models.Template.PS_CARBON_MAP_NUM}]C[C@H](C[*:{models.Template.WC_MAP_NUM_1}])[CH1:{models.Template.OLIGOMERIZATION_MAP_NUM}]=O',
    'template_pictet_spangler_kekule': None,
    'pyrroloindoline_kekule': f'[*:{models.Template.WC_MAP_NUM_1}]/C=C/[CH3:{models.Template.EAS_MAP_NUM}]',
    'aldehyde_cyclization_kekule': f'[O:{models.Template.PS_OXYGEN_MAP_NUM}]=[CH1:{models.Template.PS_CARBON_MAP_NUM}]C[C@H](C[*:{models.Template.WC_MAP_NUM_1}])[CH1:{models.Template.OLIGOMERIZATION_MAP_NUM}]=O'}
TEST_TEMPLATE_3 = {'binary': Chem.MolFromSmiles(
    'C#CCCC[C@@](Cc1cc(/C=C/COC(OC(C)(C)C)=O)ccc1)(C=O)CC(ON2C(CCC2=O)=O)=O').ToBinary(),
    'kekule': 'C#CCCC[C@@](Cc1cc(/C=C/COC(OC(C)(C)C)=O)ccc1)(C=O)CC(ON2C(CCC2=O)=O)=O',
    'oligomerization_kekule': f'C#CCCC[C@](C=O)(CC1=CC=CC(/C=C/C)=C1)C[CH:{models.Template.OLIGOMERIZATION_MAP_NUM}]=O',
    'friedel_crafts_kekule': f'[*:{models.Template.WC_MAP_NUM_1}]/C=C/[CH3:{models.Template.EAS_MAP_NUM}]',
    'tsuji_trost_kekule': f'[*:{models.Template.WC_MAP_NUM_1}]/C=C/[CH3:{models.Template.EAS_MAP_NUM}]',
    'pictet_spangler_kekule': f'C#CCCC[C@](C[*:{models.Template.WC_MAP_NUM_1}])(C[CH1:{models.Template.OLIGOMERIZATION_MAP_NUM}]=O)[CH1:{models.Template.PS_CARBON_MAP_NUM}]=[O:{models.Template.PS_OXYGEN_MAP_NUM}]',
    'template_pictet_spangler_kekule': f'C#CCCC[C@@]([CH1:{models.Template.PS_CARBON_MAP_NUM}]=[O:{models.Template.PS_OXYGEN_MAP_NUM}])(CC(=O)[NH1:{models.Template.TEMPLATE_PS_NITROGEN_MAP_NUM}][*:{models.Template.WC_MAP_NUM_1}])CC1=[CH1:{models.Template.EAS_MAP_NUM}]C=CC([*:{models.Template.WC_MAP_NUM_2}])=C1',
    'pyrroloindoline_kekule': f'[*:{models.Template.WC_MAP_NUM_1}]/C=C/[CH3:{models.Template.EAS_MAP_NUM}]',
    'aldehyde_cyclization_kekule': None}
TEST_TEMPLATES = [TEST_TEMPLATE_1, TEST_TEMPLATE_2, TEST_TEMPLATE_3]

TEST_SIDECHAIN_1 = {'binary': Chem.MolFromSmiles('CC1=CC=C(O)C=C1').ToBinary(
), 'kekule': 'CC1=CC=C(O)C=C1', 'connection': 'C', 'shared_id': 'a'}
TEST_SIDECHAIN_2 = {'binary': Chem.MolFromSmiles('CC1=COC2=NC(=O)[NH]C=C12').ToBinary(
), 'kekule': 'CC1=COC2=NC(=O)[NH]C=C12', 'connection': 'C', 'shared_id': 's'}
TEST_SIDECHAIN_3 = {'binary': Chem.MolFromSmiles('CC1=CC(=O)C2=C([NH]1)SC=C2').ToBinary(
), 'kekule': 'CC1=CC(=O)C2=C([NH]1)SC=C2', 'connection': 'C', 'shared_id': 'q'}
TEST_SIDECHAIN_4 = {'binary': Chem.MolFromSmiles('CC1=CC=CC2=N[NH]C(=O)N12').ToBinary(
), 'kekule': 'CC1=CC=CC2=N[NH]C(=O)N12', 'connection': 'C', 'shared_id': 'af'}
TEST_SIDECHAINS = [TEST_SIDECHAIN_1, TEST_SIDECHAIN_2, TEST_SIDECHAIN_3, TEST_SIDECHAIN_4]


TEST_MONOMER_1 = {'binary': Chem.MolFromSmiles('O=C(O)[C@@H]1C[C@H](OC2=CC=CC=C2)CN1').ToBinary(
), 'kekule': 'O=C(O)[C@@H]1C[C@H](OC2=CC=CC=C2)CN1', 'index': None, 'backbone': {'_id': 'alpha', 'kekule': 'NCC(=O)O'},
    'sidechain': None, 'connection': None, 'imported': True}
TEST_MONOMER_2 = {'binary': Chem.MolFromSmiles('COC1=CC=C2C(O[C@@H]3CN[C@H](C(=O)O)C3)=CC(C3=CC=CC=C3)=NC2=C1').ToBinary(
), 'kekule': 'COC1=CC=C2C(O[C@@H]3CN[C@H](C(=O)O)C3)=CC(C3=CC=CC=C3)=NC2=C1', 'index': None, 'backbone': {'_id': 'alpha', 'kekule': 'NCC(=O)O'},
    'sidechain': None, 'connection': None, 'imported': True}
TEST_MONOMER_3 = {'binary': Chem.MolFromSmiles('NC(CC1=CC=CC2=N[NH]C(=O)N12)C(=O)O').ToBinary(
), 'kekule': 'NC(CC1=CC=CC2=N[NH]C(=O)N12)C(=O)O', 'index': None, 'backbone': {'_id': 'alpha', 'kekule': 'NCC(=O)O'},
    'sidechain': 'af', 'connection': 'C', 'imported': False}
TEST_MONOMER_4 = {'binary': Chem.MolFromSmiles('NC(CC(=O)O)CC1=CC=CC2=N[NH]C(=O)N12').ToBinary(
), 'kekule': 'NC(CC(=O)O)CC1=CC=CC2=N[NH]C(=O)N12', 'index': None, 'backbone': {'_id': 'beta2', 'kekule': 'NCCC(=O)O'},
    'sidechain': 'af', 'connection': 'C', 'imported': False}
TEST_MONOMER_5 = {'binary': Chem.MolFromSmiles('NCC(CC1=CC=CC2=N[NH]C(=O)N12)C(=O)O').ToBinary(
), 'kekule': 'NCC(CC1=CC=CC2=N[NH]C(=O)N12)C(=O)O', 'index': None, 'backbone': {'_id': 'beta3', 'kekule': 'NCCC(=O)O'},
    'sidechain': 'af', 'connection': 'C', 'imported': False}
TEST_MONOMER_6 = {'binary': Chem.MolFromSmiles('NCCCC[C@H](N)C(=O)O').ToBinary(), 'kekule': 'NCCCC[C@H](N)C(=O)O',
                  'index': None, 'backbone': {'_id': 'alpha', 'kekule': 'NCC(=O)O'}, 'sidechain': None,
                  'connection': None, 'imported': True}
TEST_MONOMERS = [TEST_MONOMER_1, TEST_MONOMER_2, TEST_MONOMER_3, TEST_MONOMER_4, TEST_MONOMER_5, TEST_MONOMER_6]

TEST_PEPTIDE_1 = {'binary': Chem.MolFromSmiles(
    'NCC(=O)NC(CC(=O)NC(CC(=O)NCC(CC1=CC2=C(OC=C2)S1)C(=O)NC(CC1=CC=C2C=CC=CC=C21)C(=O)O)CC1=CC=CO1)CC1=C2C=CSC2=NS1').ToBinary(),
    'kekule': 'NCC(=O)NC(CC(=O)NC(CC(=O)NCC(CC1=CC2=C(OC=C2)S1)C(=O)NC(CC1=CC=C2C=CC=CC=C21)C(=O)O)CC1=CC=CO1)CC1=C2C=CSC2=NS1',
    'has_cap': False,
    'length': 3,
    'monomers': [
        {'_id': '12898afefgfad', 'sidechain': 'adwi8', 'proline': False},
        {'_id': 'awfseg4', 'sidechain': '3gdfbv', 'proline': False},
        {'_id': 'asfg43', 'sidechain': 'dws2', 'proline': True}
]}
TEST_PEPTIDE_2 = {'binary': Chem.MolFromSmiles('NC(CC1=CC=CC2=N[NH]C(=O)N12)C(=O)NC(CC(=O)NCC(CC1=CC=CC2=N[NH]C(=O)N12)C(=O)O)CC1=CC=CC2=N[NH]C(=O)N12').ToBinary(
), 'kekule': 'NC(CC1=CC=CC2=N[NH]C(=O)N12)C(=O)NC(CC(=O)NCC(CC1=CC=CC2=N[NH]C(=O)N12)C(=O)O)CC1=CC=CC2=N[NH]C(=O)N12',
    'has_cap': False,
    'length': 3,
    'monomers': [
    {'_id': '98asfh', 'sidechain': 'af', 'proline': False},
    {'_id': 'adjha82', 'sidechain': 'af', 'proline': False},
    {'_id': 'admaiof7', 'sidechain': 'af', 'proline': False}
]}
TEST_PEPTIDE_3 = {'binary': Chem.MolFromSmiles('COC1=CC=C2C(O[C@H]3C[C@@H](C(=O)O)N(C(=O)C(CC4=CC=CC5=N[NH]C(=O)N45)NC(=O)[C@@H]4C[C@H](OC5=CC(C6=CC=CC=C6)=NC6=CC(OC)=CC=C56)CN4C(=O)[C@@H]4C[C@H](OC5=CC=CC=C5)CN4)C3)=CC(C3=CC=CC=C3)=NC2=C1').ToBinary(
), 'kekule': 'COC1=CC=C2C(O[C@H]3C[C@@H](C(=O)O)N(C(=O)C(CC4=CC=CC5=N[NH]C(=O)N45)NC(=O)[C@@H]4C[C@H](OC5=CC(C6=CC=CC=C6)=NC6=CC(OC)=CC=C56)CN4C(=O)[C@@H]4C[C@H](OC5=CC=CC=C5)CN4)C3)=CC(C3=CC=CC=C3)=NC2=C1',
    'has_cap': False,
    'length': 4,
    'monomers': [
    {'_id': 'ad98fh', 'sidechain': None, 'proline': True},
    {'_id': 'sdwd89cvh', 'sidechain': None, 'proline': True},
    {'_id': '98asfh', 'sidechain': 'af', 'proline': False},
    {'_id': 'sdwd89cvh', 'sidechain': None, 'proline': True}
]}
TEST_PEPTIDE_4_CAP = {'binary': Chem.MolFromSmiles('COC1=CC=C2C(O[C@H]3C[C@@H](C(=O)NCCC4=CC=CC5=N[NH]C(=O)N45)N(C(=O)C(CC4=CC=CC5=N[NH]C(=O)N45)NC(=O)[C@@H]4C[C@H](OC5=CC(C6=CC=CC=C6)=NC6=CC(OC)=CC=C56)CN4C(=O)[C@@H]4C[C@H](OC5=CC=CC=C5)CN4)C3)=CC(C3=CC=CC=C3)=NC2=C1').ToBinary(
), 'kekule': 'COC1=CC=C2C(O[C@H]3C[C@@H](C(=O)NCCC4=CC=CC5=N[NH]C(=O)N45)N(C(=O)C(CC4=CC=CC5=N[NH]C(=O)N45)NC(=O)[C@@H]4C[C@H](OC5=CC(C6=CC=CC=C6)=NC6=CC(OC)=CC=C56)CN4C(=O)[C@@H]4C[C@H](OC5=CC=CC=C5)CN4)C3)=CC(C3=CC=CC=C3)=NC2=C1',
    'has_cap': True,
    'length': 4,
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
    'length': 5,
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
    'length': 5,
    'monomers': [
    {'_id': 'ad98fh', 'sidechain': None, 'proline': True},
    {'_id': 'sdwd89cvh', 'sidechain': None, 'proline': True},
    {'_id': '98asfh', 'sidechain': 'af', 'proline': False},
    {'_id': 'sdwd89cvh', 'sidechain': None, 'proline': True},
    {'_id': '12fakwd', 'sidechain': None, 'proline': False}
]}
TEST_PEPTIDES_LEN_3 = [TEST_PEPTIDE_1, TEST_PEPTIDE_2]
TEST_PEPTIDES_LEN_4 = [TEST_PEPTIDE_3, TEST_PEPTIDE_4_CAP]
TEST_PEPTIDES_LEN_5 = [TEST_PEPTIDE_4, TEST_PEPTIDE_5]

TEST_PEPTIDE_1_WITH_ID = deepcopy(TEST_PEPTIDE_1)
TEST_PEPTIDE_5_WITH_ID = deepcopy(TEST_PEPTIDE_5)
TEST_PEPTIDE_3_WITH_ID = deepcopy(TEST_PEPTIDE_3)
TEST_PEPTIDE_1_WITH_ID['_id'] = 'aefoi249'
TEST_PEPTIDE_5_WITH_ID['_id'] = 'fasef00'
TEST_PEPTIDE_3_WITH_ID['_id'] = str(uuid.uuid4())
TEST_PEPTIDE_1_WITH_ID.pop('binary')
TEST_PEPTIDE_5_WITH_ID.pop('binary')
TEST_PEPTIDE_3_WITH_ID.pop('binary')
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
TEST_TEMPLATE_PEPTIDE_8 = {'binary': Chem.MolFromSmiles('C#CC[C@@]1(C(=O)N2[C@@H]3CC[C@@]2(C(=O)NC(CCN2C=C4C=CC=CC4=C2)CC(=O)O)C3)CCCN1C(=O)C(CC1=CC=CC2=CC=CN12)NC(=O)C(CC(C)C)NC(=O)[C@H](CC=O)CC1=CC(/C=C/C)=CC=C1F').ToBinary(),
                           'kekule': 'C#CC[C@@]1(C(=O)N2[C@@H]3CC[C@@]2(C(=O)NC(CCN2C=C4C=CC=CC4=C2)CC(=O)O)C3)CCCN1C(=O)C(CC1=CC=CC2=CC=CN12)NC(=O)C(CC(C)C)NC(=O)[C@H](CC=O)CC1=CC(/C=C/C)=CC=C1F',
                           'template': str(uuid.uuid4()),
                           'peptide': TEST_PEPTIDE_5_WITH_ID}
TEST_TEMPLATE_PEPTIDE_10 = {'binary': Chem.MolFromSmiles('COC1=CC=C2C(O[C@H]3C[C@@H](C(=O)O)N(C(=O)C(CC4=CC=CC5=N[NH]C(=O)N45)NC(=O)[C@@H]4C[C@H](OC5=CC(C6=CC=CC=C6)=NC6=CC(OC)=CC=C56)CN4C(=O)[C@@H]4C[C@H](OC5=CC=CC=C5)CN4)C3)=CC(C3=CC=CC=C3)=NC2=C1').ToBinary(),
                            'kekule': 'COC1=CC=C2C(O[C@H]3C[C@@H](C(=O)O)N(C(=O)C(CC4=CC=CC5=N[NH]C(=O)N45)NC(=O)[C@@H]4C[C@H](OC5=CC(C6=CC=CC=C6)=NC6=CC(OC)=CC=C56)CN4C(=O)[C@@H]4C[C@H](OC5=CC=CC=C5)CN4)C3)=CC(C3=CC=CC=C3)=NC2=C1',
                            'template': str(uuid.uuid4()),
                            'peptide': TEST_PEPTIDE_3_WITH_ID}
TEST_TEMPLATE_PEPTIDES_LEN_3 = [TEST_TEMPLATE_PEPTIDE_1]
TEST_TEMPLATE_PEPTIDES_LEN_4 = [TEST_TEMPLATE_PEPTIDE_10]
TEST_TEMPLATE_PEPTIDES_LEN_5 = [TEST_TEMPLATE_PEPTIDE_2, TEST_TEMPLATE_PEPTIDE_3, TEST_TEMPLATE_PEPTIDE_4,
                                TEST_TEMPLATE_PEPTIDE_5, TEST_TEMPLATE_PEPTIDE_6, TEST_TEMPLATE_PEPTIDE_7, TEST_TEMPLATE_PEPTIDE_8]

TEST_MACROCYCLE_1 = {'binary': Chem.MolFromSmiles('C#CC[C@@]12CCCN1C(=O)C(CC1=CC=CC3=CC=CN13)NC(=O)C(CC(C)C)N1C=C[C@H](CC3=C(F)C=CC(=C3)/C=C/CC3=C4C=CC=CC4=CN3CCC(CC(=O)O)NC(=O)[C@]34CC[C@H](C3)N4C2=O)C1=O').ToBinary(),
                     'kekule': 'C#CC[C@@]12CCCN1C(=O)C(CC1=CC=CC3=CC=CN13)NC(=O)C(CC(C)C)N1C=C[C@H](CC3=C(F)C=CC(=C3)/C=C/CC3=C4C=CC=CC4=CN3CCC(CC(=O)O)NC(=O)[C@]34CC[C@H](C3)N4C2=O)C1=O',
                     'has_cap': False,
                     'modifications': '',
                     'template': str(uuid.uuid4()),
                     'template_peptide': str(uuid.uuid4()),
                     'reactions': [{'_id': str(uuid.uuid4()), 'type': rxns.TemplatePictetSpangler.TYPE},
                                   {'_id': str(uuid.uuid4()), 'type': rxns.FriedelCrafts.TYPE}]}


TEST_MACROCYCLE_2 = {'binary': Chem.MolFromSmiles('C#CC[C@@]12CCCN1C(=O)[C@H](CC1=CC=CC3=CC=CN13)N(C)C(=O)[C@@H](CC(C)C)N1C=C[C@H](CC3=C(F)C=CC(=C3)/C=C/CC3=C4C=CC=CC4=CN3CC[C@@H](CC(=O)O)NC(=O)[C@]34CC[C@H](C3)N4C2=O)C1=O').ToBinary(),
                     'kekule': 'C#CC[C@@]12CCCN1C(=O)[C@H](CC1=CC=CC3=CC=CN13)N(C)C(=O)[C@@H](CC(C)C)N1C=C[C@H](CC3=C(F)C=CC(=C3)/C=C/CC3=C4C=CC=CC4=CN3CC[C@@H](CC(=O)O)NC(=O)[C@]34CC[C@H](C3)N4C2=O)C1=O',
                     'has_cap': False,
                     'modifications': '',
                     'template': str(uuid.uuid4()),
                     'template_peptide': str(uuid.uuid4()),
                     'reactions': [{'_id': str(uuid.uuid4()), 'type': rxns.AldehydeCyclization.TYPE},
                                   {'_id': str(uuid.uuid4()), 'type': rxns.FriedelCrafts.TYPE}]}

TEST_REACTION_1 = {'type': 'tsuji_trost', 'binary': AllChem.ReactionFromSmarts('(c1c([*:3])ccc([OH:4])c1.C(=C/[*:50])\\[CH3:2])>>C(=C/[*:50])\\[CH2:2][O:4]c1ccc([*:3])cc1').ToBinary(),
                   'smarts': '(c1c([*:3])ccc([OH:4])c1.C(=C/[*:50])\\[CH3:2])>>C(=C/[*:50])\\[CH2:2][O:4]c1ccc([*:3])cc1',
                   'template': 'temp1',
                   'reacting_mol': {'_id': str(uuid.uuid4()), 'kekule': 'CC1=CC=C(O)C=C'},
                   'rxn_atom_idx': 5}

TEST_REACTION_2 = {'type': rxns.AldehydeCyclization.TYPE, 'binary': AllChem.ReactionFromSmarts('(O=[C:1]([C@H](C[CH:7]=O)C[*:50])[NH:6][*:5])>>O=[C:1]1[C@@H](C[*:50])C=[CH:7][N:6]1[*:5]').ToBinary(),
                   'smarts': '(O=[C:1]([C@H](C[CH:7]=O)C[*:50])[NH:6][*:5])>>O=[C:1]1[C@@H](C[*:50])C=[CH:7][N:6]1[*:5]',
                   'template': TEST_TEMPLATE_2,
                   'reacting_mol': None,
                   'rxn_atom_idx': None}

TEST_REACTION_3 = {'type': rxns.FriedelCrafts.TYPE, 'binary': AllChem.ReactionFromSmarts('(c1cc2cn([*:3])[cH:4]c2cc1.C(=C/[*:50])\\[CH3:2])>>C(=C/[*:50])\\[CH2:2][c:4]1c2ccccc2cn1[*:3]').ToBinary(),
                   'smarts': '(c1cc2cn([*:3])[cH:4]c2cc1.C(=C/[*:50])\\[CH3:2])>>C(=C/[*:50])\\[CH2:2][c:4]1c2ccccc2cn1[*:3]',
                   'template': TEST_TEMPLATE_2,
                   'reacting_mol': {'_id': str(uuid.uuid4()), 'kekule': 'CN1C=C2C=CC=CC2=C1'},
                   'rxn_atom_idx': 2}
TEST_REACTIONS = [TEST_REACTION_1, TEST_REACTION_2, TEST_REACTION_3]

TEST_REGIOSQM_PREDICTION_1 = {'predictions': [3, 6],
                              'reacting_mol': 'CC1=CC=C(O)C=C1',
                              'solvent': 'nitromethane',
                              'cutoff': '3.0'}
TEST_REGIOSQM_PREDICTION_2 = {'predictions': [2, 3, 4],
                              'reacting_mol': 'CC1=CC=C[NH]1',
                              'solvent': 'nitromethane',
                              'cutoff': '3.0'}
TEST_REGIOSQM_PREDICTIONS = [TEST_REGIOSQM_PREDICTION_1, TEST_REGIOSQM_PREDICTION_2]

TEST_PKA_PREDICTION_1 = {'predictions': {'5': 9.3},
                         'reacting_mol': 'CC1=CC=C(O)C=C1',
                         'solvent': 'water'}
TEST_PKA_PREDICTION_2 = {'predictions': {'3': 11.3, '8': 10.3},
                         'reacting_mol': 'CC1=C[NH]C2=C1C(=O)[NH]C=N2',
                         'solvent': 'water'}
TEST_PKA_PREDICTION_3 = {'predictions': {'5': 16.3},
                         'reacting_mol': 'CC1=CC=C[NH]1',
                         'solvent': 'water'}
TEST_PKA_PREDICTIONS = [TEST_PKA_PREDICTION_1, TEST_PKA_PREDICTION_2, TEST_PKA_PREDICTION_3]

FC_RESULT_SMARTS_1 = [f'(Oc1[cH:{rxns.FriedelCrafts.NUCLEOPHILE_EAS_MAP_NUM}]cc([*:{rxns.FriedelCrafts.NUCLEOPHILE_WC_MAP_NUM}])cc1.C(=C/[*:{models.Template.WC_MAP_NUM_1}])\\[CH3:{models.Template.EAS_MAP_NUM}])>>Oc1ccc([*:{rxns.FriedelCrafts.NUCLEOPHILE_WC_MAP_NUM}])c[c:{rxns.FriedelCrafts.NUCLEOPHILE_EAS_MAP_NUM}]1[CH2:{models.Template.EAS_MAP_NUM}]/C=C/[*:{models.Template.WC_MAP_NUM_1}]']
FC_RESULT_SMARTS_2 = [f'(O=[C:{rxns.FriedelCrafts.BACKBONE_CARBOXYL_MAP_NUM}]([C@@H]1C[C@H](Oc2cc[cH:{rxns.FriedelCrafts.NUCLEOPHILE_EAS_MAP_NUM}]cc2)CN1[*:{rxns.FriedelCrafts.N_TERM_WILDCARD_MAP_NUM}])[*:{rxns.FriedelCrafts.C_TERM_WILDCARD_MAP_NUM}].C(=C/[*:{models.Template.WC_MAP_NUM_1}])\\[CH3:{models.Template.EAS_MAP_NUM}])>>O=[C:{rxns.FriedelCrafts.BACKBONE_CARBOXYL_MAP_NUM}]([C@@H]1C[C@H](Oc2cc[c:{rxns.FriedelCrafts.NUCLEOPHILE_EAS_MAP_NUM}]([CH2:{models.Template.EAS_MAP_NUM}]/C=C/[*:{models.Template.WC_MAP_NUM_1}])cc2)CN1[*:{rxns.FriedelCrafts.N_TERM_WILDCARD_MAP_NUM}])[*:{rxns.FriedelCrafts.C_TERM_WILDCARD_MAP_NUM}]']
TT_RESULT_SMARTS_1 = [f'(c1c([OH:{rxns.TsujiTrost.NUCLEOPHILE_EAS_MAP_NUM}])ccc([*:{rxns.TsujiTrost.NUCLEOPHILE_WC_MAP_NUM}])c1.C(=C/[*:{models.Template.WC_MAP_NUM_1}])\\[CH3:{models.Template.EAS_MAP_NUM}])>>C(=C/[*:{models.Template.WC_MAP_NUM_1}])\\[CH2:{models.Template.EAS_MAP_NUM}][O:{rxns.TsujiTrost.NUCLEOPHILE_EAS_MAP_NUM}]c1ccc([*:{rxns.TsujiTrost.NUCLEOPHILE_WC_MAP_NUM}])cc1']
PS_RESULT_SMARTS_1 = [f'(O=[C:{rxns.PictetSpangler.TEMPLATE_OLIGOMERIZATION_MAP_NUM}]([C@@H](C[*:{models.Template.WC_MAP_NUM_1}])C[CH:{rxns.PictetSpangler.TEMPLATE_ALDEHYDE_C_MAP_NUM}]=O)[NH:{rxns.PictetSpangler.BACKBONE_NITROGEN_MAP_NUM}]C(Cc1[nH]cc[cH:{rxns.PictetSpangler.NUCLEOPHILE_EAS_MAP_NUM}]1)[*:{rxns.PictetSpangler.C_TERM_WILDCARD_MAP_NUM}])>>O=[C:{rxns.PictetSpangler.TEMPLATE_OLIGOMERIZATION_MAP_NUM}]1[C@@H](C[*:{models.Template.WC_MAP_NUM_1}])C[C:{rxns.PictetSpangler.TEMPLATE_ALDEHYDE_C_MAP_NUM}]2[c:{rxns.PictetSpangler.NUCLEOPHILE_EAS_MAP_NUM}]3cc[nH]c3CC([*:{rxns.PictetSpangler.C_TERM_WILDCARD_MAP_NUM}])[N:{rxns.PictetSpangler.BACKBONE_NITROGEN_MAP_NUM}]12']
PS_RESULT_SMARTS_2 = [f'(C#CCCC[C@@](C[C:{rxns.PictetSpangler.TEMPLATE_OLIGOMERIZATION_MAP_NUM}](=O)[NH:{rxns.PictetSpangler.BACKBONE_NITROGEN_MAP_NUM}]C(Cc1[nH]cc[cH:{rxns.PictetSpangler.NUCLEOPHILE_EAS_MAP_NUM}]1)[*:{rxns.PictetSpangler.C_TERM_WILDCARD_MAP_NUM}])(C[*:{models.Template.WC_MAP_NUM_1}])[CH:{rxns.PictetSpangler.TEMPLATE_ALDEHYDE_C_MAP_NUM}]=O)>>C#CCCC[C@]1(C[*:{models.Template.WC_MAP_NUM_1}])C[C:{rxns.PictetSpangler.TEMPLATE_OLIGOMERIZATION_MAP_NUM}](=O)[N:{rxns.PictetSpangler.BACKBONE_NITROGEN_MAP_NUM}]2C([*:{rxns.PictetSpangler.C_TERM_WILDCARD_MAP_NUM}])Cc3[nH]cc[c:{rxns.PictetSpangler.NUCLEOPHILE_EAS_MAP_NUM}]3[C:{rxns.PictetSpangler.TEMPLATE_ALDEHYDE_C_MAP_NUM}]12']
PYR_RESULT_SMARTS_1 = [f'(c1ccc2c(c1)[nH][cH:{rxns.Pyrroloindoline.ADJ_CARBON_MAP_NUM}][c:{rxns.Pyrroloindoline.NUCLEOPHILE_EAS_MAP_NUM}]2CC([NH:{rxns.Pyrroloindoline.BACKBONE_NITROGEN_MAP_NUM}][*:{rxns.Pyrroloindoline.N_TERM_WILDCARD_MAP_NUM}])[*:{rxns.Pyrroloindoline.C_TERM_WILDCARD_MAP_NUM}].C(=C/[*:{models.Template.WC_MAP_NUM_1}])\\[CH3:{rxns.Pyrroloindoline.TEMPLATE_EAS_MAP_NUM}])>>C(=C/[*:{models.Template.WC_MAP_NUM_1}])\\[CH2:{rxns.Pyrroloindoline.TEMPLATE_EAS_MAP_NUM}][C@:{rxns.Pyrroloindoline.NUCLEOPHILE_EAS_MAP_NUM}]12CC([*:{rxns.Pyrroloindoline.C_TERM_WILDCARD_MAP_NUM}])[N:{rxns.Pyrroloindoline.BACKBONE_NITROGEN_MAP_NUM}]([*:{rxns.Pyrroloindoline.N_TERM_WILDCARD_MAP_NUM}])[C@@H:{rxns.Pyrroloindoline.ADJ_CARBON_MAP_NUM}]1Nc1ccccc12',
                       f'(c1ccc2c(c1)[nH][cH:{rxns.Pyrroloindoline.ADJ_CARBON_MAP_NUM}][c:{rxns.Pyrroloindoline.NUCLEOPHILE_EAS_MAP_NUM}]2CC([NH:{rxns.Pyrroloindoline.BACKBONE_NITROGEN_MAP_NUM}][*:{rxns.Pyrroloindoline.N_TERM_WILDCARD_MAP_NUM}])[*:{rxns.Pyrroloindoline.C_TERM_WILDCARD_MAP_NUM}].C(=C/[*:{models.Template.WC_MAP_NUM_1}])\\[CH3:{rxns.Pyrroloindoline.TEMPLATE_EAS_MAP_NUM}])>>C(=C/[*:{models.Template.WC_MAP_NUM_1}])\\[CH2:{rxns.Pyrroloindoline.TEMPLATE_EAS_MAP_NUM}][C@@:{rxns.Pyrroloindoline.NUCLEOPHILE_EAS_MAP_NUM}]12CC([*:{rxns.Pyrroloindoline.C_TERM_WILDCARD_MAP_NUM}])[N:{rxns.Pyrroloindoline.BACKBONE_NITROGEN_MAP_NUM}]([*:{rxns.Pyrroloindoline.N_TERM_WILDCARD_MAP_NUM}])[C@H:{rxns.Pyrroloindoline.ADJ_CARBON_MAP_NUM}]1Nc1ccccc12']
TPS_RESULTS_SMARTS_1 = [f'(C#CCCC[C@@](CC(=O)[NH:{rxns.TemplatePictetSpangler.NITROGEN_MAP_NUM}][*:{models.Template.WC_MAP_NUM_1}])(Cc1cc([*:{models.Template.WC_MAP_NUM_2}])cc[cH:{rxns.TemplatePictetSpangler.EAS_MAP_NUM}]1)[CH:{rxns.TemplatePictetSpangler.ALDEHYDE_C_MAP_NUM}]=O)>>C#CCCC[C@]12CC(=O)[N:{rxns.TemplatePictetSpangler.NITROGEN_MAP_NUM}]([*:{models.Template.WC_MAP_NUM_1}])[C@H:{rxns.TemplatePictetSpangler.ALDEHYDE_C_MAP_NUM}]1[c:{rxns.TemplatePictetSpangler.EAS_MAP_NUM}]1ccc([*:{models.Template.WC_MAP_NUM_2}])cc1C2']
ALD_RESULT_SMARTS_1 = [f'(O=[C:{rxns.AldehydeCyclization.TEMPLATE_OLIGOMERIZATION_MAP_NUM}]([C@@H](C[*:{models.Template.WC_MAP_NUM_1}])C[CH:{rxns.AldehydeCyclization.ALDEHYDE_C_MAP_NUM}]=O)[NH:{rxns.AldehydeCyclization.BACKBONE_NITROGEN_MAP_NUM}][*:{rxns.AldehydeCyclization.BACKBONE_CARBOXYL_MAP_NUM}])>>O=[C:{rxns.AldehydeCyclization.TEMPLATE_OLIGOMERIZATION_MAP_NUM}]1[C@@H](C[*:{models.Template.WC_MAP_NUM_1}])C=[CH:{rxns.AldehydeCyclization.ALDEHYDE_C_MAP_NUM}][N:{rxns.AldehydeCyclization.BACKBONE_NITROGEN_MAP_NUM}]1[*:{rxns.AldehydeCyclization.BACKBONE_CARBOXYL_MAP_NUM}]']

TEST_PEPTIDE_PLAN_LEN_3 = models.PeptidePlan.PeptidePlanData([(1, 2, 3), (3, 4, 5)], [(6, 7, 8, 9)], 3)
TEST_PEPTIDE_PLAN_LEN_4 = models.PeptidePlan.PeptidePlanData(
    [(1, 23, 4, 5), (1, 3, 5, 5)], [(12, 4, 5, 64, 3), (1, 3, 5, 43, 2)], 4)
TEST_PEPTIDE_PLAN_LEN_5 = models.PeptidePlan.PeptidePlanData([(1, 3, 4, 5, 3), (123, 493, 123, 45, 24)], [], 5)
