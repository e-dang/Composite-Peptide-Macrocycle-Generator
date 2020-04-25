from rdkit import Chem


class Sidechain:
    TYPE = 'sidechain'

    def __init__(self, binary, kekule, connection, shared_id, _id=None):
        self._id = _id
        self.binary = binary
        self.kekule = kekule
        self.connection = connection
        self.shared_id = shared_id

    @classmethod
    def from_mol(cls, mol, connection, shared_id):
        Chem.Kekulize(mol)
        return Sidechain(mol.ToBinary(), Chem.MolToSmiles(mol, kekuleSmiles=True), connection, shared_id)

    @classmethod
    def from_dict(cls, data, _id=None):
        return Sidechain(data['binary'], data['kekule'], data['connection'], data['shared_id'], _id=_id)

    @property
    def mol(self):
        return Chem.Mol(self.binary)

    def to_dict(self):

        data = self.__dict__
        data.pop('_id')

        return data
