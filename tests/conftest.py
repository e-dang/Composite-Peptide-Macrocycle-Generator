import os
import pytest
import cpmg.hdf5 as hdf5
import cpmg.config as config
from cpmg.initializer import CPMGInitializer
import cpmg.models as models
from data.mols import *


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
def initialized_repository():
    initializer = CPMGInitializer()
    initializer.initialize()
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
def backbone_dicts(dict_sort_key):
    return sorted(TEST_BACKBONES, key=lambda x: x['mapped_kekule'])


@pytest.fixture()
def template_dicts(dict_sort_key):
    return sorted(TEST_TEMPLATES, key=dict_sort_key)


@pytest.fixture()
def sidechain_dicts(dict_sort_key):
    return sorted(TEST_SIDECHAINS, key=dict_sort_key)


@pytest.fixture()
def monomer_dicts(dict_sort_key):
    return sorted(TEST_MONOMERS, key=dict_sort_key)


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
def reaction_dicts(dict_sort_key):
    return sorted(TEST_REACTIONS, key=dict_sort_key)


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
def backbone_mols(backbone_dicts, models_from_dict):
    return models_from_dict(models.Backbone, *backbone_dicts)


@pytest.fixture()
def template_mols(template_dicts, models_from_dict):
    return models_from_dict(models.Template, *template_dicts)


@pytest.fixture()
def sidechain_mols(sidechain_dicts, models_from_dict):
    return models_from_dict(models.Sidechain, *sidechain_dicts)


@pytest.fixture()
def monomer_mols(monomer_dicts, models_from_dict):
    return models_from_dict(models.Monomer, *monomer_dicts)


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
def template_peptide_mols(template_peptide_dicts, models_from_dict):
    return models_from_dict(models.TemplatePeptide, *template_peptide_dicts)


@pytest.fixture()
def reactions(reaction_dicts, models_from_dict):
    return models_from_dict(models.Reaction, *reaction_dicts)


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
