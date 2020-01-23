import os
from abc import ABC, abstractmethod
from itertools import chain

from rdkit import Chem

import macrocycles.config as config
import macrocycles.project_io as project_io
import macrocycles.utils as utils


class IProxy(ABC):

    @abstractmethod
    def __init__(self):
        pass

    def __getitem__(self, key):

        if not self.flag:  # pylint: disable=E0203
            self.data = self.load_data()  # pylint: disable=E1103
            self.flag = True

        return self.data[key]

    @abstractmethod
    def load_data(self):
        pass


class RegioSQMProxy(IProxy):

    def __init__(self):
        self.flag = False
        self.data = {}
        self.loader = project_io.get_hashed_regiosqm_predictions

    def load_data(self):
        return self.loader()


class pKaProxy(IProxy):

    def __init__(self):
        self.flag = False
        self.data = {}
        self.loader = project_io.get_hashed_pka_predictions

    def load_data(self):
        return self.loader()


class SideChainCapProxy(IProxy):

    METHYL = Chem.MolFromSmarts('[CH3]')
    BACKBONE = Chem.MolFromSmarts('[NH2:3][CH3:2]')

    def __init__(self):
        self.flag = False
        self.data = []
        self.loader = project_io.StructureDataFileIO(os.path.join(
            config.DATA_DIR, 'chemdraw', 'conservative_sidechains.sdf'))

    def __len__(self):
        if not self.flag:
            self.data = self.load_data()
            self.flag = True

        return len(self.data)

    def load_data(self):

        return self.create_sidechain_caps(self.loader.load())

    def create_sidechain_caps(self, sidechains):

        capped_sidechains = []
        map_nums = (1, 2)
        for sidechain in sidechains:
            for atom_idx in chain.from_iterable(sidechain.GetSubstructMatches(SideChainCapProxy.METHYL)):
                atom = sidechain.GetAtomWithIdx(atom_idx)
                if atom.GetSymbol() == 'C' and atom.GetIsotope() == 0:
                    atom.SetAtomMapNum(1)
                    break

            capped_sidechains.append(utils.connect_mols(sidechain, SideChainCapProxy.BACKBONE, map_nums=map_nums))

        return capped_sidechains
