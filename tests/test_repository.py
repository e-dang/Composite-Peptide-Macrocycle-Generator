import pytest

import cpmg.repository as repo
from data.mols import *


@pytest.mark.parametrize('impl,expected_result', [(repo.HDF5, repo.hdf5.HDF5Repository)])
def test_repository_impl_from_string(impl, expected_result):

    assert repo.repository_impl_from_string(impl) is expected_result.instance()


def test_repository_impl_from_string_fail():
    with pytest.raises(ValueError):
        repo.repository_impl_from_string('dne')


def test_backbone_repository(backbone_mols, model_sort_key):
    backbone_repo = repo.create_backbone_repository()
    ids = backbone_repo.save(backbone_mols)
    data = list(backbone_repo.load(ids))

    data.sort(key=model_sort_key)

    assert data == backbone_mols


def test_backbone_repository_fail():
    backbone_repo = repo.create_backbone_repository()

    _ = backbone_repo.save(['dne'])

    assert len(backbone_repo.failed_instances) == 1


def test_connection_repository(connection_mols, model_sort_key):
    connection_repo = repo.create_connection_repository()
    ids = connection_repo.save(connection_mols)
    data = list(connection_repo.load(ids))
    data.sort(key=model_sort_key)

    assert data == connection_mols


def test_connection_repository_fail():
    connection_repo = repo.create_connection_repository()

    _ = connection_repo.save(['dne'])

    assert len(connection_repo.failed_instances) == 1


def test_template_repository(template_mols, model_sort_key):
    template_repo = repo.create_template_repository()
    ids = template_repo.save(template_mols)
    data = list(template_repo.load(ids))
    data.sort(key=model_sort_key)

    assert data == template_mols


def test_template_repository_fail():
    template_repo = repo.create_template_repository()

    _ = template_repo.save(['dne'])

    assert len(template_repo.failed_instances) == 1


def test_sidechain_repository(sidechain_mols, model_sort_key):
    sc_repo = repo.create_sidechain_repository()
    ids = sc_repo.save(sidechain_mols)
    data = list(sc_repo.load(ids))
    data.sort(key=model_sort_key)

    assert data == sidechain_mols


def test_sidechain_repository_fail():
    sc_repo = repo.create_sidechain_repository()

    _ = sc_repo.save(['dne'])

    assert len(sc_repo.failed_instances) == 1


def test_monomer_repository(monomer_mols, model_sort_key):
    monomer_repo = repo.create_monomer_repository()
    ids = monomer_repo.save(monomer_mols)
    data = list(monomer_repo.load(ids))
    data.sort(key=model_sort_key)

    for i, monomer in enumerate(monomer_mols):
        monomer.index = i

    assert data == monomer_mols


def test_monomer_repository_fail():
    monomer_repo = repo.create_monomer_repository()

    _ = monomer_repo.save(['dne'])

    assert len(monomer_repo.failed_instances) == 1


def test_peptide_repository(peptide_mols, model_sort_key):
    peptide_repo = repo.create_peptide_repository()
    ids = peptide_repo.save(peptide_mols)
    data = list(peptide_repo.load(ids))
    data.sort(key=model_sort_key)

    assert data == peptide_mols


def test_peptide_repository_fail():
    peptide_repo = repo.create_peptide_repository()

    _ = peptide_repo.save(['dne'])

    assert len(peptide_repo.failed_instances) == 1


def test_template_peptide_repository(template_peptide_mols, model_sort_key):
    template_peptide_repo = repo.create_template_peptide_repository()
    ids = template_peptide_repo.save(template_peptide_mols)
    data = list(template_peptide_repo.load(ids))
    data.sort(key=model_sort_key)

    assert data == template_peptide_mols


def test_template_peptide_repository_fail():
    template_peptide_repo = repo.create_template_peptide_repository()

    _ = template_peptide_repo.save(['dne'])

    assert len(template_peptide_repo.failed_instances) == 1


# def test_macrocycle_repository(hdf5_repository):
#     macrocycle_repo = repo.create_macrocycle_repository()
#     ids = macrocycle_repo.save(list(map(models.Macrocycle.from_dict, [TEST_MACROCYCLE_1])))
#     data = models_to_dict(macrocycle_repo.load(ids))

#     assert(len(data) == 1)
#     assert(TEST_MACROCYCLE_1 == data[0])


def test_reaction_repository(reactions, model_sort_key):
    reaction_repo = repo.create_reaction_repository()
    ids = reaction_repo.save(reactions)
    data = list(reaction_repo.load(ids))
    data.sort(key=model_sort_key)

    assert data == reactions


def test_reaction_repository_fail():
    reaction_repo = repo.create_reaction_repository()

    _ = reaction_repo.save(['dne'])

    assert len(reaction_repo.failed_instances) == 1


def test_regiosqm_repository(regiosqm_predictions, model_sort_key):
    regiosqm_repo = repo.create_regiosqm_repository()
    ids = regiosqm_repo.save(regiosqm_predictions)
    data = list(regiosqm_repo.load(ids))
    data.sort(key=model_sort_key)

    assert data == regiosqm_predictions


def test_regiosqm_repository_fail():
    regiosqm_repo = repo.create_regiosqm_repository()

    _ = regiosqm_repo.save(['dne'])

    assert len(regiosqm_repo.failed_instances) == 1


def test_pka_prediction_repository(pka_predictions, model_sort_key):
    pka_repo = repo.create_pka_repository()
    ids = pka_repo.save(pka_predictions)
    data = list(pka_repo.load(ids))
    data.sort(key=model_sort_key)

    assert data == pka_predictions


def test_pka_repository_fail():
    pka_repo = repo.create_pka_repository()

    _ = pka_repo.save(['dne'])

    assert len(pka_repo.failed_instances) == 1


def test_peptide_plan_repository(peptide_plans):
    peptide_plan_repo = repo.create_peptide_plan_repository()
    ids = peptide_plan_repo.save(peptide_plans)
    ids.set_peptide_length(peptide_plans.reg_length)

    data = peptide_plan_repo.load(ids)

    reg_combos = set(tuple(tup[1]) for tup in data.reg_combinations)
    cap_combos = set(tuple(tup[1]) for tup in data.cap_combinations)

    assert reg_combos == peptide_plans.reg_combinations
    assert cap_combos == peptide_plans.cap_combinations
    assert data.reg_length == peptide_plans.reg_length
    assert data.cap_length == peptide_plans.cap_length


def test_get_all_repository_strings():
    repo_strings = set(repo.get_all_repository_strings())

    assert repo_strings == {repo.BackboneRepository.STRING,
                            repo.ConnectionRepository.STRING,
                            repo.TemplateRepository.STRING,
                            repo.SidechainRepository.STRING,
                            repo.MonomerRepository.STRING,
                            repo.PeptideRepository.STRING,
                            repo.TemplatePeptideRepository.STRING,
                            repo.MacrocycleRepository.STRING,
                            repo.ReactionRepository.STRING,
                            repo.RegioSQMRepository.STRING,
                            repo.pKaRepository.STRING,
                            repo.PeptidePlanRepository.STRING,
                            repo.CPMGRepository.STRING}


@pytest.mark.parametrize('repository', [
    (repo.BackboneRepository),
    (repo.ConnectionRepository),
    (repo.TemplateRepository),
    (repo.SidechainRepository),
    (repo.MonomerRepository),
    (repo.PeptideRepository),
    (repo.TemplatePeptideRepository),
    (repo.MacrocycleRepository),
    (repo.ReactionRepository),
    (repo.RegioSQMRepository),
    (repo.pKaRepository),
    (repo.PeptidePlanRepository),
    (repo.CPMGRepository)
])
def test_create_repository_from_string(repository):
    produced_repository = repo.create_repository_from_string(repository.STRING)

    assert isinstance(produced_repository, repository)
