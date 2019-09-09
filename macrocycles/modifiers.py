"""
Written by Eric Dang.
github: https://github.com/e-dang
email: edang830@gmail.com
"""

from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from logging import INFO
from itertools import chain
from copy import deepcopy
import macrocycles.config as config
import macrocycles.utils as utils
from macrocycles.exceptions import CTermNotFound
import multiprocessing
from functools import partial
from collections import namedtuple

LOGGER = utils.create_logger(name=__name__, level=INFO)

NewModMacrocycle = namedtuple('NewModMacrocycle', 'mod_macrocycles parent_macrocycle')


class MacrocycleModifier(utils.Base):

    _defaults = config.DEFAULTS['MacrocycleModifier']
    AMIDE_N_MAP_NUM = 1
    METHYL_MAP_NUM = 2

    def __init__(self, logger=LOGGER, make_db_connection=True):
        """
        Initializer.

        Args:
            logger (Logger)
        """

        # I/O
        super().__init__(logger, make_db_connection)

        # data
        self.macrocycles = []

    def save_data(self):

        params = self._defaults['outputs']
        return self.to_mongo(params['col_macrocycles'])

    def load_data(self):
        try:
            params = self._defaults['inputs']
            self.macrocycles = self.from_mongo(params['col_macrocycles'], {'type': 'macrocycle'})
        except Exception:
            raise
        else:
            return True

        return False

    def generate(self, modifications=['cterm_to_amide', 'methylate']):
        func_dict = {'cterm_to_amide': MacrocycleModifier.cterm_to_amide,
                     'methylate': MacrocycleModifier.methylate}
        modifications = [func_dict[modification] for modification in modifications]

        try:
            func = partial(MacrocycleModifier.modify_macrocycle, modifications=modifications)
            with multiprocessing.Pool() as pool:
                results = pool.map_async(func, self.macrocycles)
                for mod_macrocycles, parent_macrocycle in results.get():
                    self.accumulate_data(mod_macrocycles, parent_macrocycle)
        except Exception:
            raise
        else:
            return True

        return False

    def generate_serial(self, modifications=['cterm_to_amide', 'methylate']):

        func_dict = {'cterm_to_amide': MacrocycleModifier.cterm_to_amide,
                     'methylate': MacrocycleModifier.methylate}
        modifications = [func_dict[modification] for modification in modifications]

        try:
            for macrocycle in self.macrocycles:
                mod_macrocycles, _ = MacrocycleModifier.modify_macrocycle(macrocycle, modifications)
                self.accumulate_data(mod_macrocycles, macrocycle)
        except Exception:
            raise
        else:
            return True

        return False

    def generate_from_ids(self):
        pass

    def accumulate_data(self, mod_macrocycles, parent_macrocycle):

        for i, (_, (macrocycle, modifications, kekule)) in enumerate(mod_macrocycles.items()):
            doc = deepcopy(parent_macrocycle)
            _id = ''
            for modification in modifications:
                if modification == 'amide':
                    _id += 'a'
                else:
                    _id += 'm' + str(i)

            doc['_id'] = doc['_id'] + _id
            doc['binary'] = macrocycle.ToBinary()
            doc['kekule'] = kekule
            doc['modifications'] = modifications
            self.result_data.append(doc)

        return True

    @staticmethod
    def modify_macrocycle(macrocycle, modifications):

        mol = Chem.Mol(macrocycle['binary'])
        results = {Chem.MolToSmiles(mol): [mol, [], '']}

        for modification in modifications:
            results.update(modification(deepcopy(results)))

        del results[Chem.MolToSmiles(mol)]

        return NewModMacrocycle(results, macrocycle)

    @staticmethod
    def cterm_to_amide(macrocycles):

        results = {}

        # match c-term carboxyl group
        for _, (macrocycle, modifications, _) in macrocycles.items():
            for atom in macrocycle.GetSubstructMatch(Chem.MolFromSmarts('OC(=O)[$(CN*),$(CCN*)]')):
                atom = macrocycle.GetAtomWithIdx(atom)
                if atom.GetSymbol() == 'O' and atom.GetBonds()[0].GetBondType() == Chem.BondType.SINGLE:
                    atom.SetAtomicNum(7)  # change to nitrogen
                    break
            else:
                raise CTermNotFound
            Chem.SanitizeMol(macrocycle)
            smiles = Chem.MolToSmiles(macrocycle)
            Chem.Kekulize(macrocycle)
            modifications.append('amide')
            results[smiles] = [macrocycle, modifications, Chem.MolToSmiles(macrocycle, kekuleSmiles=True)]

        return results

    @staticmethod
    def methylate(macrocycles):

        rxn = AllChem.ReactionFromSmarts('[NH1,nH1:1]>>[N,n:1]C')

        # match all amide nitrogens with one hydrogen
        results = {}
        for _, (macrocycle, modifications, _) in macrocycles.items():
            methylation_count = 0
            reactants = [macrocycle]
            while reactants:
                methylation_count += 1
                inner_results = {}
                for reactant in reactants:
                    for prod in chain.from_iterable(rxn.RunReactants((reactant,))):
                        Chem.SanitizeMol(prod)
                        smiles = Chem.MolToSmiles(prod)
                        Chem.Kekulize(prod)
                        modification = deepcopy(modifications)
                        modification.append({'methylation': methylation_count})
                        inner_results[smiles] = [prod, modification, Chem.MolToSmiles(prod, kekuleSmiles=True)]
                results.update(inner_results)
                reactants = [mol for mol, _, _ in inner_results.values()]

        return results
