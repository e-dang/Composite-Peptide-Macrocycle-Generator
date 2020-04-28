import new_architecture.repository.repository as repo
from new_architecture.io_formats import text_save
import macrocycles.config as config


class RegioSQMExporter:
    def __init__(self):
        self.sidechain_repo = repo.create_sidechain_repository()
        self.monomer_repo = repo.create_monomer_repository()

    def export_regiosqm_smiles_file(self, filepath=None):
        filepath = config.REGIOSQM_SMILES_FILEPATH if filepath is None else filepath

        data = list(filter(lambda x: x.connection != 'C', self.sidechain_repo.load()))
        data.extend(list(filter(lambda x: x.required, self.monomer_repo.load())))
        data = [str(i) + ' ' + mol.kekule + '\n' for i, mol in enumerate(data)]

        text_save(data, filepath)
