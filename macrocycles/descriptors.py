from abc import ABC
from rdkit import Chem
from rdkit.Chem import AllChem

from functools import partial
# class AbstractDescriptor(ABC):

#     def __init__(self, loader, saver, descriptor):
#         self.loader = loader
#         self.saver = saver
#         self.descriptor = descriptor
#         self.data = []

#     def run(self, **kwargs):
#         for doc in self.loader.load():
#             mol = Chem.Mol(doc['binary'])
#             self.data.append(self.descriptor(mol, **kwargs))

#         self.saver.save(self.data)


# class MolWeightDescriptor(AbstractDescriptor):

#     def __init__(self, loader, saver):
#         super().__init__(loader, saver, AllChem.CalcExactMolWt)


# class TPSADescriptor(AbstractDescriptor):

#     def __init__(self, loader, saver):
#         super().__init__(loader, saver, AllChem.CalcTPSA)

#     def run(self):
#         super().run(includeSandP=True)


# class RotatableBondsDescriptor(AbstractDescriptor):

#     def __init__(self, loader, saver):
#         super().__init__(loader, saver, AllChem.CalcNumRotatableBonds)


class MolWeightDescriptor:
    NAME = 'MW'

    def __call__(self, mol):
        return AllChem.CalcExactMolWt(mol)


class TPSADescriptor:
    NAME = 'TPSA'

    def __call__(self, mol):
        return AllChem.CalcTPSA(mol, includeSandP=True)


class RotatableBondsDescriptor:
    NAME = 'RB'

    def __call__(self, mol):
        return AllChem.CalcNumRotatableBonds(mol)
