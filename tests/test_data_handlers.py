import uuid

import pytest
import mock

import cpmg.data_handlers as handlers
import cpmg.models as models
from data.mols import *

REPO_SPEC = ['load', 'save']


def extract_ids(data):
    return sorted([doc._id for doc in data])


@pytest.fixture
def connections():
    return list(map(lambda x: models.Connection.from_dict(x[0], _id=x[1]), [(TEST_CONNECTION_1, str(uuid.uuid4())),
                                                                            (TEST_CONNECTION_2, str(uuid.uuid4()))]))


@pytest.fixture
def backbones():
    return list(map(lambda x: models.Backbone.from_dict(x[0], _id=x[1]), [(TEST_BACKBONE_1, str(uuid.uuid4())),
                                                                          (TEST_BACKBONE_2, str(uuid.uuid4())),
                                                                          (TEST_BACKBONE_3, str(uuid.uuid4()))]))


@pytest.fixture
def templates():
    return list(map(lambda x: models.Template.from_dict(x[0], _id=x[1]), [(TEST_TEMPLATE_1, str(uuid.uuid4())),
                                                                          (TEST_TEMPLATE_2, str(uuid.uuid4())),
                                                                          (TEST_TEMPLATE_3, str(uuid.uuid4()))]))


@pytest.fixture
def sidechains():
    return list(map(lambda x: models.Sidechain.from_dict(x[0], _id=x[1]), [(TEST_SIDECHAIN_1, str(uuid.uuid4())),
                                                                           (TEST_SIDECHAIN_2, str(uuid.uuid4()))]))


@pytest.fixture
def monomers():
    return list(map(lambda x: models.Monomer.from_dict(x[0], _id=x[1]), [(TEST_MONOMER_1, str(uuid.uuid4())),
                                                                         (TEST_MONOMER_4, str(uuid.uuid4()))]))


@pytest.fixture
def peptides():
    return list(map(lambda x: models.Peptide.from_dict(x[0], _id=x[1]), [(TEST_PEPTIDE_1, str(uuid.uuid4())),
                                                                         (TEST_PEPTIDE_2, str(uuid.uuid4()))]))


@pytest.fixture
def template_peptides():
    return list(map(lambda x: models.TemplatePeptide.from_dict(x[0], _id=x[1]), [(TEST_TEMPLATE_PEPTIDE_1, str(uuid.uuid4())),
                                                                                 (TEST_TEMPLATE_PEPTIDE_2, str(uuid.uuid4()))]))


@pytest.fixture
def macrocycles():
    return list(map(lambda x: models.Macrocycle.from_dict(x[0], _id=x[1]), [(TEST_MACROCYCLE_1, str(uuid.uuid4())),
                                                                            (TEST_MACROCYCLE_2, str(uuid.uuid4()))]))


@pytest.fixture
def reactions():
    return list(map(lambda x: models.Reaction.from_dict(x[0], _id=x[1]), [(TEST_REACTION_1, str(uuid.uuid4())),
                                                                          (TEST_REACTION_2, str(uuid.uuid4()))]))


@pytest.fixture
def connection_repo(connections):
    # module = request.str(getattr(module, 'handlers').__name__)
    # print(str(getattr(module, 'handlers').__name__))
    mock_connection_repo = mock.MagicMock(spec=REPO_SPEC)
    mock_connection_repo.load.return_value = connections
    mock_connection_repo.save.return_value = extract_ids(connections)
    with mock.patch('cpmg.repository.create_connection_repository', return_value=mock_connection_repo):
        yield mock_connection_repo


@pytest.fixture
def backbone_repo(backbones):
    mock_backbone_repo = mock.MagicMock(spec=REPO_SPEC)
    mock_backbone_repo.load.return_value = backbones
    mock_backbone_repo.save.return_value = extract_ids(backbones)
    with mock.patch('cpmg.repository.create_backbone_repository', return_value=mock_backbone_repo):
        yield mock_backbone_repo


@pytest.fixture
def sidechain_repo(sidechains):
    mock_sc_repo = mock.MagicMock(spec=REPO_SPEC)
    mock_sc_repo.load.return_value = sidechains
    mock_sc_repo.save.return_value = extract_ids(sidechains)
    with mock.patch('cpmg.repository.create_sidechain_repository', return_value=mock_sc_repo):
        yield mock_sc_repo


@pytest.fixture
def template_repo(templates):
    mock_template_repo = mock.MagicMock(spec=REPO_SPEC)
    mock_template_repo.load.return_value = templates
    mock_template_repo.save.return_value = extract_ids(templates)
    with mock.patch('cpmg.repository.create_template_repository', return_value=mock_template_repo):
        yield mock_template_repo


@pytest.fixture
def monomer_repo(monomers):
    mock_monomer_repo = mock.MagicMock(spec=REPO_SPEC)
    mock_monomer_repo.load.return_value = monomers
    mock_monomer_repo.save.return_value = extract_ids(monomers)
    with mock.patch('cpmg.repository.create_monomer_repository', return_value=mock_monomer_repo):
        yield mock_monomer_repo


@pytest.fixture
def peptide_repo(peptides):
    mock_peptide_repo = mock.MagicMock(spec=REPO_SPEC)
    mock_peptide_repo.load.return_value = peptides
    mock_peptide_repo.save.return_value = extract_ids(peptides)
    with mock.patch('cpmg.repository.create_peptide_repository', return_value=mock_peptide_repo):
        yield mock_peptide_repo


@pytest.fixture
def template_peptide_repo(template_peptides):
    mock_template_peptide_repo = mock.MagicMock(spec=REPO_SPEC)
    mock_template_peptide_repo.load.return_value = template_peptides
    mock_template_peptide_repo.save.return_value = extract_ids(template_peptides)
    with mock.patch('cpmg.repository.create_template_peptide_repository', return_value=mock_template_peptide_repo):
        yield mock_template_peptide_repo


@pytest.fixture
def macrocycle_repo(macrocycles):
    mock_macrocycle_repo = mock.MagicMock(spec=REPO_SPEC)
    mock_macrocycle_repo.load.return_value = macrocycles
    mock_macrocycle_repo.save.return_value = extract_ids(macrocycles)
    with mock.patch('cpmg.repository.create_macrocycle_repository', return_value=mock_macrocycle_repo):
        yield mock_macrocycle_repo


@pytest.fixture
def reaction_repo(reactions):
    mock_reaction_repo = mock.MagicMock(spec=REPO_SPEC)
    mock_reaction_repo.load.return_value = reactions
    mock_reaction_repo.save.return_value = extract_ids(reactions)
    with mock.patch('cpmg.repository.create_reaction_repository', return_value=mock_reaction_repo):
        yield mock_reaction_repo


def test_sidechain_data_handler_load(sidechain_repo, connection_repo):
    handler = handlers.SidechainDataHandler()

    data = list(handler.load())

    sidechain_repo.load.assert_called_once()
    connection_repo.load.assert_called_once()
    assert(len(data) == 2)
    for sidechain, connections in data:
        assert(isinstance(sidechain, models.Sidechain))
        assert(isinstance(connections, list))
        assert(len(connections) == 1)
        assert(connections[0].to_dict() == TEST_CONNECTION_2)


def test_sidechain_data_handler_save(sidechains, sidechain_repo, connection_repo):
    handler = handlers.SidechainDataHandler()

    ids = handler.save(sidechains)
    ids.sort()

    sidechain_repo.save.assert_called_once()
    assert(ids == extract_ids(sidechains))


def test_monomer_data_handler_load(sidechain_repo, backbone_repo):
    handler = handlers.MonomerDataHandler()

    data = list(handler.load())

    sidechain_repo.load.assert_called_once()
    backbone_repo.load.assert_called_once()
    assert(len(data) == 2)
    for sidechain, backbones in data:
        assert(isinstance(sidechain, models.Sidechain))
        assert(isinstance(backbones, list))
        assert(len(backbones) == 3)
        assert(all(map(lambda x: isinstance(x, models.Backbone), backbones)))


def test_monomer_data_handler_save(monomers, monomer_repo):
    handler = handlers.MonomerDataHandler()

    ids = handler.save(monomers)

    monomer_repo.save.assert_called_once()
    assert(ids == extract_ids(monomers))


def test_template_peptide_data_handler_load(peptide_repo, template_repo):
    handler = handlers.TemplatePeptideDataHandler()

    data = list(handler.load())

    peptide_repo.load.assert_called_once()
    template_repo.load.assert_called_once()
    assert(len(data) == 2)
    for peptide, templates in data:
        assert(isinstance(peptide, models.Peptide))
        assert(isinstance(templates, list))
        assert(len(templates) == 3)
        assert(all(map(lambda x: isinstance(x, models.Template), templates)))


def test_template_peptide_data_handler_save(template_peptides, template_peptide_repo):
    handler = handlers.TemplatePeptideDataHandler()

    ids = handler.save(template_peptides)

    template_peptide_repo.save.assert_called_once()
    assert(ids == extract_ids(template_peptides))


# def test_macrocycle_data_handler_load(template_peptide_repo, reaction_repo):
#     handler = handlers.MacrocycleDataHandler()

#     data = list(handler.load())

#     template_peptide_repo.load.assert_called_once()
#     reaction_repo.load.assert_called_once()
#     assert(len(data) == 2)
#     assert()

def test_intermolecular_reaction_data_handler_load(sidechain_repo, monomer_repo):
    handler = handlers.InterMolecularReactionDataHandler()

    data = list(handler.load())

    sidechain_repo.load.assert_called_once()
    monomer_repo.load.assert_called_once()
    assert(len(data) == 4)
    for i, doc in enumerate(data):
        if i < 2:
            assert(isinstance(doc, models.Sidechain))
        else:
            assert(isinstance(doc, models.Monomer))


def test_intremolecular_reaction_data_handler_save(reactions, reaction_repo):
    handler = handlers.InterMolecularReactionDataHandler()

    ids = handler.save(reactions)

    reaction_repo.save.assert_called_once()
    assert(ids == extract_ids(reactions))


def test_intramolecular_reaction_data_handler_load(template_repo):
    handler = handlers.IntraMolecularReactionDataHandler()

    data = list(handler.load())

    template_repo.load.assert_called_once()
    assert(len(data) == 3)
    for template in data:
        assert(isinstance(template, models.Template))


def test_intramolecular_reaction_data_handler_save(reactions, reaction_repo):
    handler = handlers.IntraMolecularReactionDataHandler()

    ids = handler.save(reactions)

    reaction_repo.save.assert_called_once()
    assert(ids == extract_ids(reactions))
