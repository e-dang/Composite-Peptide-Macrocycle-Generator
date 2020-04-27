from rdkit import Chem
from rdkit.Chem import AllChem
from macrocycles.exceptions import InvalidMolecule
import new_architecture.utils as utils
import macrocycles.exceptions as exceptions
from itertools import chain

SC_ATTACHMENT_POINT = Chem.MolFromSmarts('[CH3;!13CH3]')  # methyls marked with C13 aren't used as attachment points
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
    MAP_NUM = 1

    def __init__(self, binary, kekule, mapped_kekule, _id=None):
        super().__init__(binary, kekule, _id)
        self.mapped_kekule = mapped_kekule

    def __eq__(self, other):
        return self.kekule == other.kekule and self.mapped_kekule == other.mapped_kekule

    @classmethod
    def from_mol(cls, mol):
        Chem.SanitizeMol(mol)
        cls.validate(mol)
        Chem.Kekulize(mol)
        mapped_kekule = Chem.MolToSmiles(mol, kekuleSmiles=True)
        utils.clear_atom_map_nums(mol)
        return cls(mol.ToBinary(), Chem.MolToSmiles(mol, kekuleSmiles=True), mapped_kekule)

    @classmethod
    def from_dict(cls, data, _id=None):
        cls.validate(Chem.Mol(data['binary']))
        return cls(data['binary'], data['kekule'], data['mapped_kekule'], _id=_id)

    @staticmethod
    def validate(mol):
        try:
            _, map_nums = zip(*utils.get_atom_map_nums(mol))
            if not all(map(lambda x: x == Backbone.MAP_NUM, map_nums)):
                raise ValueError
        except ValueError:
            raise exceptions.InvalidMolecule(
                f'Backbone molecule is missing an atom map number or the atom map number is not equal to 1')

        return True

    def to_reduced_dict(self):
        data = self.__dict__
        data.pop('binary')
        data.pop('mapped_kekule')
        return data


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
    OLIGOMERIZATION_MAP_NUM = 1

    def __init__(self, binary, kekule, mapped_kekule, _id=None):
        super().__init__(binary, kekule, _id)
        self.mapped_kekule = mapped_kekule

    @classmethod
    def from_mol(cls, mol, full_kekule):
        Chem.Kekulize(mol)
        cls.validate(mol)
        return cls(mol.ToBinary(), full_kekule, Chem.MolToSmiles(mol, kekuleSmiles=True))

    @classmethod
    def from_dict(cls, data, _id=None):
        cls.validate(Chem.Mol(data['binary']))
        return cls(data['binary'], data['kekule'], data['mapped_kekule'], _id=_id)

    @staticmethod
    def validate(mol):
        try:
            _, map_nums = zip(*utils.get_atom_map_nums(mol))
            if Template.OLIGOMERIZATION_MAP_NUM not in map_nums:
                raise ValueError(f'Template molecule is missing oligomerization atom map number!')
        except ValueError as err:
            raise exceptions.InvalidMolecule(str(err))

        return True


class Sidechain(AbstractMolecule):
    MAP_NUM = 2

    def __init__(self, binary, kekule, attachment_point, connection, shared_id, _id=None):
        super().__init__(binary, kekule, _id)
        self.attachment_point = attachment_point
        self.connection = connection
        self.shared_id = shared_id

    def __eq__(self, other):
        return self.kekule == other.kekule and self.connection == other.connection and self.shared_id == other.shared_id

    @classmethod
    def from_mol(cls, mol, connection, shared_id):
        Chem.SanitizeMol(mol)
        Chem.Kekulize(mol)
        attachment_point = cls.validate(mol)
        return cls(mol.ToBinary(), Chem.MolToSmiles(mol, kekuleSmiles=True), attachment_point, connection._id, shared_id)

    @classmethod
    def from_dict(cls, data, _id=None):
        attachment_point = cls.validate(Chem.Mol(data['binary']))
        return cls(data['binary'], data['kekule'], attachment_point, data['connection'], data['shared_id'], _id=_id)

    @staticmethod
    def validate(mol):
        attachment_point = list(chain.from_iterable(mol.GetSubstructMatches(SC_ATTACHMENT_POINT)))
        if len(attachment_point) != 1:
            raise exceptions.InvalidMolecule(f'Sidechains must have exactly one attachment point')

        return attachment_point[0]

    @property
    def mapped_mol(self):
        sidechain = self.mol
        sidechain.GetAtomWithIdx(self.attachment_point).SetAtomMapNum(self.MAP_NUM)
        return sidechain

    def to_dict(self):
        data = super().to_dict()
        data.pop('attachment_point')
        return data


class Monomer(AbstractMolecule):
    def __init__(self, binary, kekule, required, backbone, sidechain, connection, proline, imported, _id=None, index=None):
        super().__init__(binary, kekule, _id)
        self.index = index
        self.required = required
        self.backbone = backbone
        self.sidechain = sidechain
        self.connection = connection
        self.proline = proline
        self.imported = imported

    def __eq__(self, other):
        return self.kekule == other.kekule and self.backbone == other.backbone and self.sidechain == other.sidechain \
            and self.required == other.required and self.proline == other.proline and self.imported == other.imported

    @classmethod
    def from_mol(cls, mol, backbone, sidechain, imported=False):
        Chem.Kekulize(mol)
        return cls(mol.ToBinary(), Chem.MolToSmiles(mol, kekuleSmiles=True), cls.is_required(mol),
                   backbone.to_reduced_dict(), sidechain.shared_id, sidechain.connection, cls.is_proline(mol), imported)

    @classmethod
    def from_dict(cls, data, _id=None):
        mol = Chem.Mol(data['binary'])
        return cls(data['binary'], data['kekule'], cls.is_required(mol), data['backbone'], data['sidechain'],
                   data['connection'], cls.is_proline(mol), data['imported'], _id=_id, index=data['index'])

    @staticmethod
    def is_required(mol):
        return bool(AllChem.CalcNumAromaticRings(mol))

    @staticmethod
    def is_proline(mol):
        return bool(AllChem.CalcNumAliphaticRings(mol) and mol.HasSubstructMatch(PROLINE_N_TERM))

    @property
    def backbone_mol(self):
        return Chem.MolFromSmiles(self.backbone['kekule'])

    def to_dict(self):
        data = super().to_dict()
        data.pop('required')
        data.pop('proline')
        return data


class Peptide(AbstractMolecule):
    def __init__(self, binary, kekule, has_cap, monomers, _id=None):
        super().__init__(binary, kekule, _id)
        self.has_cap = has_cap
        self.monomers = monomers

    def __eq__(self, other):
        return self.kekule == other.kekule and self.has_cap == other.has_cap and self.monomers == other.monomers

    @classmethod
    def from_mol(cls, mol, has_cap, monomers):
        Chem.Kekulize(mol)
        monomers = [{'_id': monomer._id, 'sidechain': monomer.sidechain,
                     'proline': monomer.proline} for monomer in monomers]
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
