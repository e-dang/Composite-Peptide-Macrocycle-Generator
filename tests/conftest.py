import os
import uuid
from copy import deepcopy

import pytest

import cpmg.config as config
import cpmg.hdf5 as hdf5
import cpmg.models as models
from cpmg.initializer import CPMGInitializer
from cpmg.exporters import RegioSQMExporter
from data.mols import *


def add_ids(data):
    data = deepcopy(data)
    for doc in data:
        doc['_id'] = str(uuid.uuid4())
    return data


def add_ids_to_peptide_plan_tuple(plan_tup):
    data = []
    for tup in plan_tup.reg_combos + plan_tup.cap_combos:
        data.append((str(uuid.uuid4()), tup))

    return data


@pytest.fixture(autouse=True)
def patched_configs(tmpdir, monkeypatch):
    monkeypatch.setattr('cpmg.config.DATA_DIR', os.path.join(config.TEST_DIR, 'data'))
    monkeypatch.setattr('cpmg.config.IMPORT_DIR', os.path.join(config.TEST_DIR, 'data', 'imports'))
    monkeypatch.setattr('cpmg.config.DATA_FORMAT', 'hdf5')
    monkeypatch.setattr('cpmg.config.HDF5_FILEPATH', os.path.join(str(tmpdir), 'test_hdf5_repo.hdf5'))
    monkeypatch.setattr('cpmg.config.REGIOSQM_SMILES_FILEPATH', os.path.join(str(tmpdir), 'regiosqm_smiles.txt'))


@pytest.fixture(autouse=True)
def hdf5_reset():
    yield None
    if hasattr(hdf5.HDF5Repository, '_instance'):
        del hdf5.HDF5Repository._instance


@pytest.fixture()
def partial_initialized_repository():
    initializer = CPMGInitializer()
    initializer.initialize_mols_only()
    yield None


@pytest.fixture()
def initialized_repository(partial_initialized_repository):
    exporter = RegioSQMExporter()
    exporter.export_regiosqm_smiles_file()
    initializer = CPMGInitializer()
    initializer.initialize_predictions_only()
    yield None


@pytest.fixture()
def dict_sort_key():
    def func(x):
        if x.get('mapped_kekule') is not None:
            return x['mapped_kekule']

        if x.get('smarts') is not None:
            return x['smarts']

        if x.get('reacting_mol') is not None:
            return x['reacting_mol']

        return x['kekule']

    return func


@pytest.fixture()
def model_sort_key():
    def func(x):
        if isinstance(x, models.Backbone):
            return x.mapped_kekule

        if isinstance(x, models.Reaction):
            return x.smarts

        if isinstance(x, models.AbstractPrediction):
            return x.reacting_mol

        return x.kekule

    return func


@pytest.fixture(autouse=True)
def models_from_dict():
    def func(model, *args):
        try:
            return [model.from_dict(doc, _id=doc['_id']) for doc in args]
        except KeyError:
            return [model.from_dict(doc) for doc in args]

    return func


@pytest.fixture(autouse=True)
def models_to_dict():
    def func(model):
        return list(map(lambda x: x.to_dict(), model))

    return func


@pytest.fixture(autouse=True)
def extract_attr_from_models():
    def func(data, attr):
        return sorted([getattr(model, attr) for model in data])

    return func


@pytest.fixture(autouse=True)
def extract_attr_from_dicts():
    def func(docs, attr):
        return sorted([doc[attr] for doc in docs])

    return func


@pytest.fixture()
def connection_dicts(dict_sort_key):
    return sorted(TEST_CONNECTIONS, key=dict_sort_key)


@pytest.fixture()
def connection_w_id_dicts(dict_sort_key):
    return sorted(add_ids(TEST_CONNECTIONS), key=dict_sort_key)


@pytest.fixture()
def backbone_dicts(dict_sort_key):
    return sorted(TEST_BACKBONES, key=lambda x: x['mapped_kekule'])


@pytest.fixture()
def backbone_w_id_dicts(dict_sort_key):
    return sorted(add_ids(TEST_BACKBONES), key=lambda x: x['mapped_kekule'])


@pytest.fixture()
def template_dicts(dict_sort_key):
    return sorted(TEST_TEMPLATES, key=dict_sort_key)


@pytest.fixture()
def template_w_id_dicts(dict_sort_key):
    return sorted(add_ids(TEST_TEMPLATES), key=dict_sort_key)


@pytest.fixture()
def sidechain_dicts(dict_sort_key):
    return sorted(TEST_SIDECHAINS, key=dict_sort_key)


@pytest.fixture()
def sidechain_w_id_dicts(dict_sort_key):
    return sorted(add_ids(TEST_SIDECHAINS), key=dict_sort_key)


@pytest.fixture()
def monomer_dicts(dict_sort_key):
    return sorted(TEST_MONOMERS, key=dict_sort_key)


@pytest.fixture()
def monomer_w_id_dicts(dict_sort_key):
    return sorted(add_ids(TEST_MONOMERS_W_IDXS), key=dict_sort_key)


@pytest.fixture()
def peptide_len_3_dicts(dict_sort_key):
    return sorted(TEST_PEPTIDES_LEN_3, key=dict_sort_key)


@pytest.fixture()
def peptide_len_4_dicts(dict_sort_key):
    return sorted(TEST_PEPTIDES_LEN_4, key=dict_sort_key)


@pytest.fixture()
def peptide_len_5_dicts(dict_sort_key):
    return sorted(TEST_PEPTIDES_LEN_5, key=dict_sort_key)


@pytest.fixture()
def peptide_dicts(dict_sort_key):
    return sorted(TEST_PEPTIDES_LEN_3 + TEST_PEPTIDES_LEN_4 + TEST_PEPTIDES_LEN_5, key=dict_sort_key)


@pytest.fixture()
def peptide_w_id_dicts(dict_sort_key):
    return sorted(add_ids(TEST_PEPTIDES_LEN_3 + TEST_PEPTIDES_LEN_4 + TEST_PEPTIDES_LEN_5), key=dict_sort_key)


@pytest.fixture()
def template_peptide_len_3_dicts(dict_sort_key):
    return sorted(TEST_TEMPLATE_PEPTIDES_LEN_3, key=dict_sort_key)


@pytest.fixture()
def template_peptide_len_4_dicts(dict_sort_key):
    return sorted(TEST_TEMPLATE_PEPTIDES_LEN_4, key=dict_sort_key)


@pytest.fixture()
def template_peptide_len_5_dicts(dict_sort_key):
    return sorted(TEST_TEMPLATE_PEPTIDES_LEN_5, key=dict_sort_key)


@pytest.fixture()
def template_peptide_dicts(dict_sort_key):
    return sorted(TEST_TEMPLATE_PEPTIDES_LEN_3 + TEST_TEMPLATE_PEPTIDES_LEN_4 + TEST_TEMPLATE_PEPTIDES_LEN_5, key=dict_sort_key)


@pytest.fixture()
def template_peptide_w_id_dicts(dict_sort_key):
    return sorted(add_ids(TEST_TEMPLATE_PEPTIDES_LEN_3 + TEST_TEMPLATE_PEPTIDES_LEN_4 + TEST_TEMPLATE_PEPTIDES_LEN_5), key=dict_sort_key)


@pytest.fixture()
def macrocycle_w_id_dicts(dict_sort_key):
    return sorted(add_ids(TEST_MACROCYCLES), key=dict_sort_key)


@pytest.fixture()
def reaction_dicts(dict_sort_key):
    return sorted(TEST_REACTIONS, key=dict_sort_key)


@pytest.fixture()
def reaction_w_id_dicts(dict_sort_key):
    return sorted(add_ids(TEST_REACTIONS), key=dict_sort_key)


@pytest.fixture()
def regiosqm_dicts(dict_sort_key):
    return sorted(TEST_REGIOSQM_PREDICTIONS, key=dict_sort_key)


@pytest.fixture()
def pka_dicts(dict_sort_key):
    return sorted(TEST_PKA_PREDICTIONS, key=dict_sort_key)


@pytest.fixture()
def peptide_plan_tuple_len_3():
    return TEST_PEPTIDE_PLAN_LEN_3


@pytest.fixture()
def peptide_plan_tuple_len_4():
    return TEST_PEPTIDE_PLAN_LEN_4


@pytest.fixture()
def peptide_plan_tuple_len_5():
    return TEST_PEPTIDE_PLAN_LEN_5


@pytest.fixture()
def connection_mols(connection_dicts, models_from_dict):
    return models_from_dict(models.Connection, *connection_dicts)


@pytest.fixture()
def connection_mols_w_ids(connection_w_id_dicts, models_from_dict):
    return models_from_dict(models.Connection, *connection_w_id_dicts)


@pytest.fixture()
def backbone_mols(backbone_dicts, models_from_dict):
    return models_from_dict(models.Backbone, *backbone_dicts)


@pytest.fixture()
def backbone_w_id_mols(backbone_w_id_dicts, models_from_dict):
    return models_from_dict(models.Backbone, *backbone_w_id_dicts)


@pytest.fixture()
def template_mols(template_dicts, models_from_dict):
    return models_from_dict(models.Template, *template_dicts)


@pytest.fixture()
def template_w_id_mols(template_w_id_dicts, models_from_dict):
    return models_from_dict(models.Template, *template_w_id_dicts)


@pytest.fixture()
def sidechain_mols(sidechain_dicts, models_from_dict):
    return models_from_dict(models.Sidechain, *sidechain_dicts)


@pytest.fixture()
def sidechain_w_id_mols(sidechain_w_id_dicts, models_from_dict):
    return models_from_dict(models.Sidechain, *sidechain_w_id_dicts)


@pytest.fixture()
def monomer_mols(monomer_dicts, models_from_dict):
    return models_from_dict(models.Monomer, *monomer_dicts)


@pytest.fixture()
def monomer_w_idx_mols(monomer_w_id_dicts, models_from_dict):
    return models_from_dict(models.Monomer, *monomer_w_id_dicts)


@pytest.fixture()
def peptide_len_3_mols(peptide_len_3_dicts, models_from_dict):
    return models_from_dict(models.Peptide, *peptide_len_3_dicts)


@pytest.fixture()
def peptide_len_4_mols(peptide_len_4_dicts, models_from_dict):
    return models_from_dict(models.Peptide, *peptide_len_4_dicts)


@pytest.fixture()
def peptide_len_5_mols(peptide_len_5_dicts, models_from_dict):
    return models_from_dict(models.Peptide, *peptide_len_5_dicts)


@pytest.fixture()
def peptide_mols(peptide_dicts, models_from_dict):
    return models_from_dict(models.Peptide, *peptide_dicts)


@pytest.fixture()
def peptide_w_id_mols(peptide_w_id_dicts, models_from_dict):
    return models_from_dict(models.Peptide, *peptide_w_id_dicts)


@pytest.fixture()
def template_peptide_mols(template_peptide_dicts, models_from_dict):
    return models_from_dict(models.TemplatePeptide, *template_peptide_dicts)


@pytest.fixture()
def template_peptide_w_id_mols(template_peptide_w_id_dicts, models_from_dict):
    return models_from_dict(models.TemplatePeptide, *template_peptide_w_id_dicts)


@pytest.fixture()
def macrocycle_w_id_mols(macrocycle_w_id_dicts, models_from_dict):
    return models_from_dict(models.Macrocycle, *macrocycle_w_id_dicts)


@pytest.fixture()
def reactions(reaction_dicts, models_from_dict):
    return models_from_dict(models.Reaction, *reaction_dicts)


@pytest.fixture()
def reactions_w_ids(reaction_w_id_dicts, models_from_dict):
    return models_from_dict(models.Reaction, *reaction_w_id_dicts)


@pytest.fixture()
def regiosqm_predictions(regiosqm_dicts, models_from_dict):
    return models_from_dict(models.RegioSQMPrediction, *regiosqm_dicts)


@pytest.fixture()
def pka_predictions(pka_dicts, models_from_dict):
    return models_from_dict(models.pKaPrediction, *pka_dicts)


@pytest.fixture(params=['peptide_plan_tuple_len_3', 'peptide_plan_tuple_len_4', 'peptide_plan_tuple_len_5'])
def peptide_plans(request):
    plan_tup = request.getfixturevalue(request.param)
    peptide_plan = models.PeptidePlan(plan_tup.length)

    for tup in plan_tup.reg_combos + plan_tup.cap_combos:
        peptide_plan.add(tup)

    return peptide_plan


@pytest.fixture(params=['peptide_plan_tuple_len_3', 'peptide_plan_tuple_len_4', 'peptide_plan_tuple_len_5'])
def peptide_plans_w_ids(request):
    plan_tup = request.getfixturevalue(request.param)
    peptide_plan = models.PeptidePlan.from_array_tuple(plan_tup.length, add_ids_to_peptide_plan_tuple(plan_tup))
    return peptide_plan
