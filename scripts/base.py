from pathlib import Path
from utils import Database


class Base():

    def __init__(self, project_dir=Path(__file__).resolve().parents[1], db_rxn=False, db_mol=False):

        self.project_dir = project_dir

        if db_rxn:
            self.db_rxn = Database(host='localhost', port=270127, db='rxn_templates')

        if db_mol:
            self.db_mol = Database(host='localhost', port=270127, db='molecules')
