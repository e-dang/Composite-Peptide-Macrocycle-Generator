import cpmg.repository as repo


def test_initializer(initialized_repository):
    assert(len(list(repo.create_backbone_repository().load())) == 3)
    assert(len(list(repo.create_connection_repository().load())) == 2)
    assert(len(list(repo.create_template_repository().load())) == 3)
    assert(len(list(repo.create_sidechain_repository().load())) == 3)
    assert(len(list(repo.create_monomer_repository().load())) == 4)
