from rdkit import Chem
from rdkit.Chem import AllChem
from macrocycles.exceptions import InvalidMolecule
import new_architecture.utils as utils
import macrocycles.exceptions as exceptions

NON_ATTACHMENT_METHYL = Chem.MolFromSmarts('[CH3;!13CH3]')  # methyls marked with C13 aren't used as attachment points
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

# class AbstractReaction:
#     def __init__(self, binary, smarts, _id=None):
#         self.binary = binary
#         self.smarts = smarts
#         self.rxn_atom_idx = rxn_atom_idx
#         self.template = template
#         self.reacting_mol = reacting_mol


class Backbone(AbstractMolecule):
    def __init__(self, binary, kekule, mapped_kekule, _id=None):
        super().__init__(binary, kekule, _id)
        self.mapped_kekule = mapped_kekule

    @classmethod
    def from_mol(cls, mol):
        cls.validate(mol)
        Chem.Kekulize(mol)
        mapped_kekule = Chem.MolToSmiles(mol, kekuleSmiles=True)
        utils.clear_atom_map_nums(mol)
        return cls(mol.ToBinary(), Chem.MolToSmiles(mol, kekuleSmiles=True), mapped_kekule)

    @classmethod
    def from_dict(cls, data, _id=None):
        cls.validate(Chem.MolFromSmiles(data['mapped_kekule']))
        return cls(data['binary'], data['kekule'], data['mapped_kekule'], _id=_id)

    @staticmethod
    def validate(mol):
        try:
            _, map_nums = zip(*utils.get_atom_map_nums(mol))
            if not all(map(lambda x: x == 1, map_nums)):
                raise ValueError
        except ValueError:
            raise exceptions.InvalidMolecule(
                f'Backbone molecule is missing an atom map number or the atom map number is not equal to 1')

        return True


class Connection(AbstractMolecule):
    def __init__(self, binary, kekule, _id=None):
        super().__init__(binary, kekule, _id)

    @classmethod
    def from_mol(cls, mol):
        Chem.Kekulize(mol)
        return cls(mol.ToBinary(), Chem.MolToSmiles(mol, kekuleSmiles=True))

    @classmethod
    def from_dict(cls, data, _id=None):
        return cls(data['binary'], data['kekule'], _id=_id)


class Template(AbstractMolecule):
    def __init__(self, binary, kekule, _id=None):
        super().__init__(binary, kekule, _id)

    @classmethod
    def from_mol(cls, mol):
        Chem.Kekulize(mol)
        return cls(mol.ToBinary(), Chem.MolToSmiles(mol, kekuleSmiles=True))

    @classmethod
    def from_dict(cls, data, _id=None):
        return cls(data['binary'], data['kekule'], _id=_id)


class Sidechain(AbstractMolecule):
    def __init__(self, binary, kekule, connection, shared_id, _id=None):
        super().__init__(binary, kekule, _id)
        self.connection = connection
        self.shared_id = shared_id

    def __eq__(self, other):
        return self.kekule == other.kekule and self.connection == other.connection and self.shared_id == other.shared_id

    @classmethod
    def from_mol(cls, mol, connection, shared_id):
        Chem.Kekulize(mol)
        return cls(mol.ToBinary(), Chem.MolToSmiles(mol, kekuleSmiles=True), connection._id, shared_id)

    @classmethod
    def from_dict(cls, data, _id=None):
        return cls(data['binary'], data['kekule'], data['connection'], data['shared_id'], _id=_id)


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
    def from_mol(cls, mol, backbone, sidechain, imported=False):
        Chem.Kekulize(mol)
        required = bool(AllChem.CalcNumAromaticRings(mol))
        is_proline = bool(AllChem.CalcNumAliphaticRings(mol) and mol.HasSubstructMatch(PROLINE_N_TERM))
        return cls(mol.ToBinary(), Chem.MolToSmiles(mol, kekuleSmiles=True), required, backbone._id,
                   sidechain.shared_id, sidechain.connection, is_proline, imported)

    @classmethod
    def from_dict(cls, data, _id=None):
        return cls(data['binary'], data['kekule'], data['required'], data['backbone'], data['sidechain'],
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
        return cls(mol.ToBinary(), Chem.MolToSmiles(mol, kekuleSmiles=True), has_cap, monomers)

    @classmethod
    def from_dict(cls, data, _id=None):
        return cls(data['binary'], data['kekule'], data['has_cap'], data['monomers'], _id=_id)


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
        return cls(mol.ToBinary(), Chem.MolToSmiles(mol, kekuleSmiles=True), template._id, peptide)

    @classmethod
    def from_dict(cls, data, _id=None):
        return cls(data['binary'], data['kekule'], data['template'], data['peptide'], _id=_id)


class Macrocycle(AbstractMolecule):
    def __init__(self, binary, kekule, modifications, has_cap, template_peptide, template, reactions, _id=None):
        super().__init__(binary, kekule, _id)
        self.modifications = modifications
        self.has_cap = has_cap
        self.template_peptide = template_peptide
        self.template = template
        self.reactions = reactions

    @classmethod
    def from_mol(cls, mol, modifications, template_peptide, reactions):
        Chem.Kekulize(mol)
        reactions = [rxn._id for rxn in reactions]
        return cls(mol.ToBinary(), Chem.MolToSmiles(mol, kekuleSmiles=True), modifications,
                   template_peptide.has_cap, template_peptide._id, template_peptide.template, reactions)

    @classmethod
    def from_dict(cls, data, _id=None):
        return cls(data['binary'], data['kekule'], data['modifications'], data['has_cap'],
                   data['template_peptide'], data['template'], data['reactions'], _id=_id)
