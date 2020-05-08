import uuid

import pytest
import mock

import cpmg.ranges as ranges
import cpmg.data_handlers as handlers
import cpmg.models as models
from data.mols import *

REPO_SPEC = ['load', 'save']


def extract_ids(data):
    return sorted([doc._id for doc in data])


@pytest.fixture
def connection_repo(connection_mols_w_ids):
    mock_connection_repo = mock.MagicMock(spec=REPO_SPEC)
    mock_connection_repo.load.return_value = connection_mols_w_ids
    mock_connection_repo.save.return_value = extract_ids(connection_mols_w_ids)
    with mock.patch('cpmg.repository.create_connection_repository', return_value=mock_connection_repo):
        yield mock_connection_repo


@pytest.fixture
def backbone_repo(backbone_w_id_mols):
    mock_backbone_repo = mock.MagicMock(spec=REPO_SPEC)
    mock_backbone_repo.load.return_value = backbone_w_id_mols
    mock_backbone_repo.save.return_value = extract_ids(backbone_w_id_mols)
    with mock.patch('cpmg.repository.create_backbone_repository', return_value=mock_backbone_repo):
        yield mock_backbone_repo


@pytest.fixture
def sidechain_repo(sidechain_w_id_mols):
    mock_sc_repo = mock.MagicMock(spec=REPO_SPEC)
    mock_sc_repo.load.return_value = sidechain_w_id_mols
    mock_sc_repo.save.return_value = extract_ids(sidechain_w_id_mols)
    with mock.patch('cpmg.repository.create_sidechain_repository', return_value=mock_sc_repo):
        yield mock_sc_repo


@pytest.fixture
def template_repo(template_w_id_mols):
    mock_template_repo = mock.MagicMock(spec=REPO_SPEC)
    mock_template_repo.load.return_value = template_w_id_mols
    mock_template_repo.save.return_value = extract_ids(template_w_id_mols)
    with mock.patch('cpmg.repository.create_template_repository', return_value=mock_template_repo):
        yield mock_template_repo


@pytest.fixture
def monomer_repo(monomer_w_idx_mols):
    mock_monomer_repo = mock.MagicMock(spec=REPO_SPEC)
    mock_monomer_repo.load.return_value = monomer_w_idx_mols
    mock_monomer_repo.save.return_value = extract_ids(monomer_w_idx_mols)
    with mock.patch('cpmg.repository.create_monomer_repository', return_value=mock_monomer_repo):
        yield mock_monomer_repo


@pytest.fixture
def peptide_repo(peptide_w_id_mols):
    mock_peptide_repo = mock.MagicMock(spec=REPO_SPEC)
    mock_peptide_repo.load.return_value = peptide_w_id_mols
    mock_peptide_repo.save.return_value = extract_ids(peptide_w_id_mols)
    with mock.patch('cpmg.repository.create_peptide_repository', return_value=mock_peptide_repo):
        yield mock_peptide_repo


@pytest.fixture
def template_peptide_repo(template_peptide_w_id_mols):
    mock_template_peptide_repo = mock.MagicMock(spec=REPO_SPEC)
    mock_template_peptide_repo.load.return_value = template_peptide_w_id_mols
    mock_template_peptide_repo.save.return_value = extract_ids(template_peptide_w_id_mols)
    with mock.patch('cpmg.repository.create_template_peptide_repository', return_value=mock_template_peptide_repo):
        yield mock_template_peptide_repo


@pytest.fixture
def macrocycle_repo(macrocycle_w_id_mols):
    mock_macrocycle_repo = mock.MagicMock(spec=REPO_SPEC)
    mock_macrocycle_repo.load.return_value = macrocycle_w_id_mols
    mock_macrocycle_repo.save.return_value = extract_ids(macrocycle_w_id_mols)
    with mock.patch('cpmg.repository.create_macrocycle_repository', return_value=mock_macrocycle_repo):
        yield mock_macrocycle_repo


@pytest.fixture
def reaction_repo(reactions_w_ids):
    mock_reaction_repo = mock.MagicMock(spec=REPO_SPEC)
    mock_reaction_repo.load.return_value = reactions_w_ids
    mock_reaction_repo.save.return_value = extract_ids(reactions_w_ids)
    with mock.patch('cpmg.repository.create_reaction_repository', return_value=mock_reaction_repo):
        yield mock_reaction_repo


@pytest.fixture
def peptide_plan_repo(peptide_plans):
    mock_peptide_plan_repo = mock.MagicMock(spec=REPO_SPEC)
    mock_peptide_plan_repo.load.return_value = peptide_plans
    mock_peptide_plan_repo.save.return_value = extract_ids(peptide_plans)
    with mock.patch('cpmg.repository.create_peptide_plan_repository', return_value=mock_peptide_plan_repo):
        yield mock_peptide_plan_repo, peptide_plans


def test_sidechain_data_handler_load(sidechain_repo, connection_repo, sidechain_w_id_mols):
    handler = handlers.SidechainDataHandler()

    data = list(handler.load())

    sidechain_repo.load.assert_called_once()
    connection_repo.load.assert_called_once()
    assert(len(data) == len(sidechain_w_id_mols))
    for sidechain, connections in data:
        assert(isinstance(sidechain, models.Sidechain))
        assert(isinstance(connections, list))
        assert(len(connections) == 1)
        assert(connections[0].to_dict() == TEST_CONNECTION_2)


def test_sidechain_data_handler_save(sidechain_w_id_mols, sidechain_repo, connection_repo):
    handler = handlers.SidechainDataHandler()

    ids = handler.save(sidechain_w_id_mols)
    ids.sort()

    sidechain_repo.save.assert_called_once()
    assert(ids == extract_ids(sidechain_w_id_mols))


def test_monomer_data_handler_load(sidechain_repo, backbone_repo, sidechain_w_id_mols, backbone_w_id_mols):
    handler = handlers.MonomerDataHandler()

    data = list(handler.load())

    sidechain_repo.load.assert_called_once()
    backbone_repo.load.assert_called_once()
    assert(len(data) == len(sidechain_w_id_mols))
    for sidechain, backbones in data:
        assert(isinstance(sidechain, models.Sidechain))
        assert(isinstance(backbones, list))
        assert(len(backbones) == len(backbone_w_id_mols))
        assert(all(map(lambda x: isinstance(x, models.Backbone), backbones)))


def test_monomer_data_handler_save(monomer_w_idx_mols, monomer_repo):
    handler = handlers.MonomerDataHandler()

    ids = handler.save(monomer_w_idx_mols)

    monomer_repo.save.assert_called_once()
    assert(ids == extract_ids(monomer_w_idx_mols))


# def test_peptide_data_handler_load(monomer_repo, peptide_plan_repo):
#     plan_repo, peptide_plan = peptide_plan_repo
#     handler = handlers.PeptideDataHandler()

#     data = list(handler.load()

#     monomer_repo.load.assert_called_once()
#     plan_repo.load.assert_called_once()
#     assert len(data) == len(peptide_plan)
#     indices = []
#     for monomer_combination in data:
#         combo_indices = []
#         for monomer in monomer_combination:
#             assert isinstance(monomer, models.Monomer)
#             combo_indices.append(monomer.index)
#         indices.append(tuple(combo_indices))


def test_template_peptide_data_handler_load(peptide_repo, template_repo, peptide_w_id_mols, template_w_id_mols):
    handler = handlers.TemplatePeptideDataHandler()

    data = list(handler.load())

    peptide_repo.load.assert_called_once()
    template_repo.load.assert_called_once()
    assert(len(data) == len(peptide_w_id_mols))
    for peptide, templates in data:
        assert(isinstance(peptide, models.Peptide))
        assert(isinstance(templates, list))
        assert(len(templates) == len(template_w_id_mols))
        assert(all(map(lambda x: isinstance(x, models.Template), templates)))


def test_template_peptide_data_handler_save(template_peptide_w_id_mols, template_peptide_repo):
    handler = handlers.TemplatePeptideDataHandler()

    ids = handler.save(template_peptide_w_id_mols)

    template_peptide_repo.save.assert_called_once()
    assert(ids == extract_ids(template_peptide_w_id_mols))


# def test_macrocycle_data_handler_load(template_peptide_repo, reaction_repo):
#     handler = handlers.MacrocycleDataHandler()

#     data = list(handler.load())

#     template_peptide_repo.load.assert_called_once()
#     reaction_repo.load.assert_called_once()
#     assert(len(data) == 2)
#     assert()

def test_intermolecular_reaction_data_handler_load(sidechain_repo, monomer_repo, sidechain_w_id_mols, monomer_w_idx_mols):
    handler = handlers.InterMolecularReactionDataHandler()

    data = list(handler.load())

    sidechain_repo.load.assert_called_once()
    monomer_repo.load.assert_called_once()
    assert len(data) == len(sidechain_w_id_mols) + len(monomer_w_idx_mols)
    for i, doc in enumerate(data):
        if i < len(sidechain_w_id_mols):
            assert(isinstance(doc, models.Sidechain))
        else:
            assert(isinstance(doc, models.Monomer))


def test_intremolecular_reaction_data_handler_save(reactions_w_ids, reaction_repo):
    handler = handlers.InterMolecularReactionDataHandler()

    ids = handler.save(reactions_w_ids)

    reaction_repo.save.assert_called_once()
    assert(ids == extract_ids(reactions_w_ids))


def test_intramolecular_reaction_data_handler_load(template_repo):
    handler = handlers.IntraMolecularReactionDataHandler()

    data = list(handler.load())

    template_repo.load.assert_called_once()
    assert(len(data) == 3)
    for template in data:
        assert(isinstance(template, models.Template))


def test_intramolecular_reaction_data_handler_save(reactions_w_ids, reaction_repo):
    handler = handlers.IntraMolecularReactionDataHandler()

    ids = handler.save(reactions_w_ids)

    reaction_repo.save.assert_called_once()
    assert(ids == extract_ids(reactions_w_ids))
