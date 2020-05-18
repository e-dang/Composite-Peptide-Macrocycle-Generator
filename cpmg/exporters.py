import cpmg.repository as repo
from cpmg.utils import save_text
import cpmg.config as config


class RegioSQMExporter:
    def __init__(self):
        self.sidechain_repo = repo.create_sidechain_repository()
        self.monomer_repo = repo.create_monomer_repository()

    def export_regiosqm_smiles_file(self, filepath=None):
        filepath = filepath or config.REGIOSQM_SMILES_FILEPATH

        data = list(filter(lambda x: x.connection == 'C', self.sidechain_repo.load()))
        data.extend(list(filter(lambda x: x.required, self.monomer_repo.load())))
        data = [str(i) + ' ' + mol.kekule + '\n' for i, mol in enumerate(data)]

        save_text(data, filepath)

        return True
