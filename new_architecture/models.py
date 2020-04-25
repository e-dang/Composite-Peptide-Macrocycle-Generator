from rdkit import Chem
from rdkit.Chem import AllChem

PROLINE_N_TERM = Chem.MolFromSmarts('[NH;R]')


class AbstractMolecule:
    def __init__(self, binary, kekule, _id=None):
        self._id = _id
        self.binary = binary
        self.kekule = kekule

    @property
    def mol(self):
        return Chem.Mol(self.binary)

    def to_dict(self):
        data = self.__dict__
        data.pop('_id')
        return data


class Backbone(AbstractMolecule):
    def __init__(self, binary, kekule, _id=None):
        super().__init__(binary, kekule, _id)

    @classmethod
    def from_mol(cls, mol):
        Chem.Kekulize(mol)
        return Connection(mol.ToBinary(), Chem.MolToSmiles(mol, kekuleSmiles=True))

    @classmethod
    def from_dict(cls, data, _id=None):
        return Connection(data['binary'], data['kekule'], _id=_id)


class Connection(AbstractMolecule):
    def __init__(self, binary, kekule, _id=None):
        super().__init__(binary, kekule, _id)

    @classmethod
    def from_mol(cls, mol):
        Chem.Kekulize(mol)
        return Connection(mol.ToBinary(), Chem.MolToSmiles(mol, kekuleSmiles=True))

    @classmethod
    def from_dict(cls, data, _id=None):
        return Connection(data['binary'], data['kekule'], _id=_id)


class Template(AbstractMolecule):
    def __init__(self, binary, kekule, _id=None):
        super().__init__(binary, kekule, _id)

    @classmethod
    def from_mol(cls, mol):
        Chem.Kekulize(mol)
        return Template(mol.ToBinary(), Chem.MolToSmiles(mol, kekuleSmiles=True))

    @classmethod
    def from_dict(cls, data, _id=None):
        return Template(data['binary'], data['kekule'], _id=_id)


class Sidechain(AbstractMolecule):
    def __init__(self, binary, kekule, connection, shared_id, _id=None):
        super().__init__(binary, kekule, _id)
        self.connection = connection
        self.shared_id = shared_id

    @classmethod
    def from_mol(cls, mol, connection, shared_id):
        Chem.Kekulize(mol)
        return Sidechain(mol.ToBinary(), Chem.MolToSmiles(mol, kekuleSmiles=True), connection, shared_id)

    @classmethod
    def from_dict(cls, data, _id=None):
        return Sidechain(data['binary'], data['kekule'], data['connection'], data['shared_id'], _id=_id)


class Monomer(AbstractMolecule):
    def __init__(self, binary, kekule, required, backbone, sidechain, connection, is_proline, imported, _id=None, index=None):
        super().__init__(binary, kekule, _id)
        self.index = index
        self.required = required
        self.backbone = backbone
        self.sidechain = sidechain
        self.connection = connection
        self.is_proline = is_proline
        self.imported = imported

    @classmethod
    def from_mol(cls, mol, backbone, sidechain):
        Chem.Kekulize(mol)
        required = bool(AllChem.CalcNumAromaticRings(mol))
        is_proline = AllChem.CalcNumAliphaticRings(mol) and mol.HasSubstructMatch(PROLINE_N_TERM)
        return Monomer(mol.ToBinary(), Chem.MolToSmiles(mol, kekuleSmiles=True), required, backbone,
                       sidechain.shared_id, sidechain.connection, is_proline, False)

    @classmethod
    def from_dict(cls, data, _id=None):
        return Monomer(data['binary'], data['kekule'], data['required'], data['backbone'], data['sidechain'],
                       data['connection'], data['is_proline'], data['imported'], _id=_id, index=data['index'])


class Peptide(AbstractMolecule):
    def __init__(self, binary, kekule, has_cap, monomers, _id=None):
        super().__init__(binary, kekule, _id)
        self.has_cap = has_cap
        self.monomers = monomers

    @classmethod
    def from_mol(cls, mol, has_cap, monomers):
        Chem.Kekulize(mol)
        monomers = [{'_id': monomer._id, 'sidechain': monomer.sidechain,
                     'is_proline': monomer.is_proline} for monomer in monomers]
        return Peptide(mol.ToBinary(), Chem.MolToSmiles(mol, kekuleSmiles=True), has_cap, monomers)

    @classmethod
    def from_dict(cls, data, _id=None):
        return Peptide(data['binary'], data['kekule'], data['has_cap'], data['monomers'], _id=_id)


class TemplatePeptide(AbstractMolecule):
    def __init__(self, binary, kekule, template, peptide, _id=None):
        super().__init__(binary, kekule, _id)
        self.template = template
        self.peptide = peptide

    @classmethod
    def from_mol(cls, mol, template, peptide):
        Chem.Kekulize(mol)
        peptide = peptide.__dict__
        peptide.pop('binary')
        return TemplatePeptide(mol.ToBinary(), Chem.MolToSmiles(mol, kekuleSmiles=True), template._id, peptide)

    @classmethod
    def from_dict(cls, data, _id=None):
        return TemplatePeptide(data['binary'], data['kekule'], data['template'], data['peptide'], _id=_id)
