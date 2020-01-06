from abc import ABC, abstractmethod

from rdkit import Chem
from rdkit.Chem import AllChem
import config
import iterators
import project_io


class IImporter(ABC):

    @abstractmethod
    def import_data(self):
        pass


class DataImporter(IImporter):

    def __init__(self, data_format=config.DATA_FORMAT):
        if data_format == 'json':
            sidechain_io = project_io.JsonSideChainIO()
            monomer_io = project_io.JsonMonomerIO()
            id_io = project_io.JsonIDIO()
            index_io = project_io.JsonIndexIO()
            regiosqm_io = project_io.JsonRegioSQMIO()
            pka_io = project_io.JsonpKaIO()
        elif data_format == 'mongo':
            sidechain_io = project_io.MongoSideChainIO()
            monomer_io = project_io.MongoMonomerIO()
            id_io = project_io.MongoIDIO()
            index_io = project_io.MongoIndexIO()
            regiosqm_io = project_io.MongoRegioSQMIO()
            pka_io = project_io.MongopKaIO()

        id_iterator = iterators.IDIterator(id_io)
        index_iterator = iterators.IndexIterator(index_io)
        self.molecule_importers = [SideChainImporter(sidechain_io, id_iterator),
                                   MonomerImporter(monomer_io, id_iterator, index_iterator)]
        self.prediction_importers = [RegioSQMImporter(regiosqm_io, sidechain_io), pKaImporter(pka_io, sidechain_io)]

    def import_data(self):
        self.import_molecules()
        self.import_predictions()

    def import_molecules(self):
        for importer in self.molecule_importers:
            importer.import_data()

    def import_predictions(self):
        for importer in self.prediction_importers:
            importer.import_data()


class SideChainImporter(IImporter):

    def __init__(self, sidechain_io, id_iterator):

        self.loader = project_io.ChemDrawSideChainIO()
        self.saver = sidechain_io
        self.id_iterator = id_iterator

    def import_data(self):
        data = []
        for sidechain in self.loader.load():
            _id = self.id_iterator.get_next()
            sidechain = Chem.MolFromSmiles(Chem.MolToSmiles(sidechain))  # canonicalize
            binary = sidechain.ToBinary()
            Chem.Kekulize(sidechain)
            data.append({'_id': _id,
                         'type': 'sidechain',
                         'binary': binary,
                         'kekule': Chem.MolToSmiles(sidechain, kekuleSmiles=True),
                         'connection': 'methyl',
                         'shared_id': _id})

        self.saver.save(data)
        self.id_iterator.save()


class MonomerImporter(IImporter):

    def __init__(self, monomer_io, id_iterator, index_iterator):

        self.loader = project_io.ChemDrawMonomerIO()
        self.saver = monomer_io
        self.id_iterator = id_iterator
        self.index_iterator = index_iterator

    def import_data(self):

        data = []
        patt = Chem.MolFromSmarts('[NH;R]')
        for monomer in self.loader.load():
            binary = monomer.ToBinary()
            monomer = Chem.MolFromSmiles(Chem.MolToSmiles(monomer))  # canonicalize
            proline = bool(AllChem.CalcNumAliphaticRings(monomer)) and monomer.HasSubstructMatch(patt)
            Chem.Kekulize(monomer)
            data.append({'_id': self.id_iterator.get_next().upper(),
                         'type': 'monomer',
                         'binary': binary,
                         'kekule': Chem.MolToSmiles(monomer, kekuleSmiles=True),
                         'index': self.index_iterator.get_next(),
                         'required': bool(AllChem.CalcNumAromaticRings(monomer)),
                         'backbone': 'alpha',
                         'sidechain': None,
                         'proline': proline})

        self.saver.save(data)
        self.index_iterator.save()


class RegioSQMImporter(IImporter):

    def __init__(self, regiosqm_io, sidechain_io, solvent='nitromethane', cutoff=3):

        self.raw_regiosqm_io = project_io.RawRegioSQMIO()
        self.regiosqm_io = regiosqm_io
        self.sidechain_io = sidechain_io
        self.solvent = solvent
        self.cutoff = cutoff
        self.hashed_sidechains = {}

    def hash_sidechains(self):

        for sidechain in filter(lambda x: x['connection'] == 'methyl', self.sidechain_io.load()):
            self.hashed_sidechains[sidechain['_id']] = sidechain['kekule']

    def import_data(self):

        data = []
        self.hash_sidechains()
        for i, row in enumerate(self.raw_regiosqm_io.load()):
            if i % 3 == 0:
                sidechain_id = row[0]
            elif i % 3 == 1:
                predictions = [int(row[j]) for j in range(0, len(row), 2)]
            else:
                self.check(sidechain_id, predictions)
                data.append({'type': 'regiosqm',
                             'predictions': predictions,
                             'sidechain': sidechain_id,
                             'solvent': self.solvent,
                             'cutoff': self.cutoff})

        self.regiosqm_io.save(data)

    def check(self, sidechain_id, predictions):

        sidechain = Chem.MolFromSmiles(self.hashed_sidechains[sidechain_id])
        Chem.Kekulize(sidechain)
        for atom_idx in predictions:
            atom = sidechain.GetAtomWithIdx(atom_idx)
            assert atom.GetSymbol() == 'C'
            assert atom.GetTotalNumHs() > 0
            assert atom.GetIsAromatic()


class pKaImporter(IImporter):

    def __init__(self, pka_io, sidechain_io, solvent='water'):

        self.raw_pka_io = project_io.RawpKaIO()
        self.pka_io = pka_io
        self.sidechain_io = sidechain_io
        self.solvent = solvent
        self.hashed_sidechains = {}

    def hash_sidechains(self):

        for sidechain in filter(lambda x: x['connection'] == 'methyl', self.sidechain_io.load()):
            self.hashed_sidechains[sidechain['_id']] = sidechain['kekule']

    def import_data(self):

        data = []
        self.hash_sidechains()
        for row in self.raw_pka_io.load():
            _id, smiles, pkas = row.split(';')
            pkas = list(map(float, pkas.strip('\n').strip(' ').split(',')))

            # atom map numbers can change the ordering of atoms in unpredictable ways, thus we must do this VERY
            # annoying round about way of mapping the predictions to the correct heteroatoms on the sidechains used to
            # create reactions
            sidechain_w_map_nums = Chem.MolFromSmiles(smiles.strip(' '))
            sidechain = Chem.MolFromSmiles(self.hashed_sidechains[_id])
            Chem.Kekulize(sidechain_w_map_nums)
            Chem.Kekulize(sidechain)

            # map atom indexes and pkas to each other in sidechain_w_map_nums using atom map numbers
            pka_map = {}
            for atom in sidechain_w_map_nums.GetAtoms():
                map_num = atom.GetAtomMapNum()
                if map_num:
                    pka_map[map_num - 1] = (atom.GetIdx(), pkas[map_num - 1])
            map_list = [pka_map[i] for i in sorted(pka_map)]

            # match the pkas to the heteroatoms on properly ordered sidechain
            for match in sidechain.GetSubstructMatches(sidechain_w_map_nums):
                heteroatom_pkas = {str(match[atom_idx]): pka for atom_idx, pka in map_list}

            # check because canonicalization of smiles strings might not be consistent
            self.check(_id, heteroatom_pkas)

            # format data
            data.append({'type': 'pka',
                         'predictions': heteroatom_pkas,
                         'sidechain': _id,
                         'solvent': self.solvent})

        self.pka_io.save(data)

    def check(self, sidechain_id, predictions):

        sidechain = Chem.MolFromSmiles(self.hashed_sidechains[sidechain_id])
        Chem.Kekulize(sidechain)
        for atom_idx, _ in predictions.items():
            atom = sidechain.GetAtomWithIdx(int(atom_idx))
            assert atom.GetSymbol() in ('N', 'O', 'S')
            assert atom.GetTotalNumHs() > 0
