from copy import deepcopy
from itertools import chain

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

import cpmg.config as config
import cpmg.exceptions as exceptions
import cpmg.utils as utils

SC_ATTACHMENT_POINT = Chem.MolFromSmarts('[CH3;!13CH3]')  # methyls marked with C13 aren't used as attachment points
METHYL = Chem.MolFromSmarts('[CH3]')
CARBOXYL = Chem.MolFromSmarts('C(=O)[OH]')
N_TERM = Chem.MolFromSmarts('[NH2]')
PROLINE_N_TERM = Chem.MolFromSmarts('[NH;R]')
ALPHA_BACKBONE = Chem.MolFromSmarts('NCC(=O)O')
TAGGED_ALPHA_BACKBONE = Chem.MolFromSmarts('N[CH2:1]C(=O)[OH]')


class AbstractMolecule:
    def __init__(self, binary, kekule, _id=None):
        self._id = _id
        self.binary = binary
        self.kekule = kekule

    def __eq__(self, other):
        return self.kekule == other.kekule

    @property
    def mol(self):
        return Chem.Mol(self.binary)

    def to_dict(self):
        data = deepcopy(self.__dict__)
        data.pop('_id')
        return data


class Backbone(AbstractMolecule):
    MAP_NUM = 1

    def __init__(self, binary, kekule, mapped_kekule, _id=None):
        super().__init__(binary, kekule, _id)
        self.mapped_kekule = mapped_kekule

    def __eq__(self, other):
        return super().__eq__(other) and self.mapped_kekule == other.mapped_kekule

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
        data = deepcopy(self.__dict__)
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
    EAS_MAP_NUM = 200
    WC_MAP_NUM_1 = 201
    WC_MAP_NUM_2 = 202
    PS_OXYGEN_MAP_NUM = 300
    PS_CARBON_MAP_NUM = 301
    TEMPLATE_PS_NITROGEN_MAP_NUM = 302

    def __init__(self, binary, kekule, oligomerization_kekule, friedel_crafts_kekule, tsuji_trost_kekule,
                 pictet_spangler_kekule, template_pictet_spangler_kekule, pyrroloindoline_kekule,
                 aldehyde_cyclization_kekule, _id=None):
        super().__init__(binary, kekule, _id)
        self.oligomerization_kekule = oligomerization_kekule
        self.friedel_crafts_kekule = friedel_crafts_kekule
        self.tsuji_trost_kekule = tsuji_trost_kekule
        self.pictet_spangler_kekule = pictet_spangler_kekule
        self.template_pictet_spangler_kekule = template_pictet_spangler_kekule
        self.pyrroloindoline_kekule = pyrroloindoline_kekule
        self.aldehyde_cyclization_kekule = aldehyde_cyclization_kekule

    @classmethod
    def from_mol(cls, mol, oligomerization_kekule, friedel_crafts_kekule, tsuji_trost_kekule, pictet_spangler_kekule,
                 template_pictet_spangler_kekule, pyrroloindoline_kekule, aldehyde_cyclization_kekule):
        cls.validate({'oligomerization_kekule': oligomerization_kekule,
                      'friedel_crafts_kekule': friedel_crafts_kekule,
                      'tsuji_trost_kekule': tsuji_trost_kekule,
                      'pictet_spangler_kekule': pictet_spangler_kekule,
                      'template_pictet_spangler_kekule': template_pictet_spangler_kekule,
                      'pyrroloindoline_kekule': pyrroloindoline_kekule,
                      'aldehyde_cyclization_kekule': aldehyde_cyclization_kekule})
        Chem.Kekulize(mol)
        return cls(mol.ToBinary(), Chem.MolToSmiles(mol, kekuleSmiles=True), oligomerization_kekule,
                   friedel_crafts_kekule, tsuji_trost_kekule, pictet_spangler_kekule, template_pictet_spangler_kekule,
                   pyrroloindoline_kekule, aldehyde_cyclization_kekule)

    @classmethod
    def from_dict(cls, data, _id=None):
        cls.validate(data)
        return cls(data['binary'], data['kekule'], data['oligomerization_kekule'], data['friedel_crafts_kekule'],
                   data['tsuji_trost_kekule'], data['pictet_spangler_kekule'], data['template_pictet_spangler_kekule'],
                   data['pyrroloindoline_kekule'], data['aldehyde_cyclization_kekule'], _id=_id)

    @staticmethod
    def validate(data):
        try:
            if data['oligomerization_kekule'] is not None:
                Template.validate_oligomerization_mol(Chem.MolFromSmiles(data['oligomerization_kekule']))
            if data['friedel_crafts_kekule'] is not None:
                Template.validate_friedel_crafts_mol(Chem.MolFromSmiles(data['friedel_crafts_kekule']))
            if data['tsuji_trost_kekule'] is not None:
                Template.validate_tsuji_trost_mol(Chem.MolFromSmiles(data['tsuji_trost_kekule']))
            if data['pictet_spangler_kekule'] is not None:
                Template.validate_pictet_spangler_mol(Chem.MolFromSmiles(data['pictet_spangler_kekule']))
            if data['template_pictet_spangler_kekule'] is not None:
                Template.validate_template_pictet_spangler_mol(
                    Chem.MolFromSmiles(data['template_pictet_spangler_kekule']))
            if data['pyrroloindoline_kekule'] is not None:
                Template.validate_pyrroloindoline_mol(Chem.MolFromSmiles(data['pyrroloindoline_kekule']))
            if data['aldehyde_cyclization_kekule'] is not None:
                Template.validate_aldehyde_cyclization_mol(Chem.MolFromSmiles(data['aldehyde_cyclization_kekule']))
        except ValueError as err:
            raise exceptions.InvalidMolecule(str(err))

        return True

    @staticmethod
    def validate_oligomerization_mol(mol):
        _, map_nums = zip(*utils.get_atom_map_nums(mol))
        if Template.OLIGOMERIZATION_MAP_NUM not in map_nums:
            raise ValueError(f'Template molecule is missing oligomerization atom map number!')

        return True

    @staticmethod
    def validate_friedel_crafts_mol(mol):
        _, map_nums = zip(*utils.get_atom_map_nums(mol))
        if Template.EAS_MAP_NUM not in map_nums or Template.WC_MAP_NUM_1 not in map_nums:
            raise ValueError(f'Template molecule is missing friedel crafts atom map numbers!')

        return True

    @staticmethod
    def validate_tsuji_trost_mol(mol):
        _, map_nums = zip(*utils.get_atom_map_nums(mol))
        if Template.EAS_MAP_NUM not in map_nums or Template.WC_MAP_NUM_1 not in map_nums:
            raise ValueError(f'Template molecule is missing tsuji trost atom map numbers!')

        return True

    @staticmethod
    def validate_pictet_spangler_mol(mol):
        _, map_nums = zip(*utils.get_atom_map_nums(mol))
        if Template.OLIGOMERIZATION_MAP_NUM not in map_nums \
                or Template.WC_MAP_NUM_1 not in map_nums \
                or Template.PS_OXYGEN_MAP_NUM not in map_nums \
                or Template.PS_CARBON_MAP_NUM not in map_nums:
            raise ValueError(f'Template molecule is missing pictet spangler atom map numbers!')

        return True

    @staticmethod
    def validate_template_pictet_spangler_mol(mol):
        _, map_nums = zip(*utils.get_atom_map_nums(mol))
        if Template.EAS_MAP_NUM not in map_nums \
                or Template.WC_MAP_NUM_1 not in map_nums \
                or Template.WC_MAP_NUM_2 not in map_nums \
                or Template.PS_OXYGEN_MAP_NUM not in map_nums \
                or Template.PS_CARBON_MAP_NUM not in map_nums \
                or Template.TEMPLATE_PS_NITROGEN_MAP_NUM not in map_nums:
            raise ValueError(f'Template molecule is missing pictet spangler atom map numbers!')

        return True

    @staticmethod
    def validate_pyrroloindoline_mol(mol):
        _, map_nums = zip(*utils.get_atom_map_nums(mol))
        if Template.EAS_MAP_NUM not in map_nums or Template.WC_MAP_NUM_1 not in map_nums:
            raise ValueError(f'Template molecule is missing pyrroloindoline atom map numbers!')

        return True

    @staticmethod
    def validate_aldehyde_cyclization_mol(mol):
        _, map_nums = zip(*utils.get_atom_map_nums(mol))
        if Template.OLIGOMERIZATION_MAP_NUM not in map_nums \
                or Template.WC_MAP_NUM_1 not in map_nums \
                or Template.PS_OXYGEN_MAP_NUM not in map_nums \
                or Template.PS_CARBON_MAP_NUM not in map_nums:
            raise ValueError(f'Template molecule is missing pictet spangler atom map numbers!')

        return True

    @property
    def oligomerization_mol(self):
        return Chem.MolFromSmiles(self.oligomerization_kekule)

    @property
    def friedel_crafts_mol(self):
        return Chem.MolFromSmiles(self.friedel_crafts_kekule) if self.friedel_crafts_kekule is not None else None

    @property
    def tsuji_trost_mol(self):
        return Chem.MolFromSmiles(self.tsuji_trost_kekule) if self.tsuji_trost_kekule is not None else None

    @property
    def pictet_spangler_mol(self):
        return Chem.MolFromSmiles(self.pictet_spangler_kekule) if self.pictet_spangler_kekule is not None else None

    @property
    def template_pictet_spangler_mol(self):
        return Chem.MolFromSmiles(self.template_pictet_spangler_kekule) if self.template_pictet_spangler_kekule is not None else None

    @property
    def pyrroloindoline_mol(self):
        return Chem.MolFromSmiles(self.pyrroloindoline_kekule) if self.pyrroloindoline_kekule is not None else None

    @property
    def aldehyde_cyclization_mol(self):
        return Chem.MolFromSmiles(self.aldehyde_cyclization_kekule) if self.aldehyde_cyclization_kekule is not None else None


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
        return cls(mol.ToBinary(), Chem.MolToSmiles(mol, kekuleSmiles=True), attachment_point, connection.kekule, shared_id)

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

    def __eq__(self, other):
        return self.kekule == other.kekule and self.template == other.template and self.peptide == other.peptide

    @classmethod
    def from_mol(cls, mol, template, peptide):
        Chem.Kekulize(mol)
        peptide = deepcopy(peptide.__dict__)
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
        mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
        Chem.Kekulize(mol)
        cls.validate(mol)
        reactions = [rxn.to_reduced_dict() for rxn in reactions]
        return cls(mol.ToBinary(), Chem.MolToSmiles(mol, kekuleSmiles=True), modifications,
                   template_peptide.peptide['has_cap'], template_peptide._id, template_peptide.template, reactions)

    @classmethod
    def from_dict(cls, data, _id=None):
        return cls(data['binary'], data['kekule'], data['modifications'], data['has_cap'],
                   data['template_peptide'], data['template'], data['reactions'], _id=_id)

    @staticmethod
    def validate(mol):
        if not [ring for ring in mol.GetRingInfo().BondRings() if len(ring) >= config.MIN_MACRO_RING_SIZE]:
            raise exceptions.InvalidMolecule(
                f'A macrocycle must have at least {config.MIN_MACRO_RING_SIZE} ring atoms!')

        return True


class Reaction:
    def __init__(self, rxn_type, binary, smarts, template, reacting_mol, rxn_atom_idx, _id=None):
        self._id = _id
        self.type = rxn_type
        self.binary = binary
        self.smarts = smarts
        self.template = template
        self.reacting_mol = reacting_mol
        self.rxn_atom_idx = rxn_atom_idx

    def __eq__(self, other):
        return self.type == other.type and self.smarts == other.smarts and self.template == other.template \
            and self.reacting_mol == other.reacting_mol and self.rxn_atom_idx == other.rxn_atom_idx

    @classmethod
    def from_mols(cls, rxn_type, smarts, template, reacting_mol, rxn_atom_idx):

        if reacting_mol is not None:
            _id = reacting_mol.shared_id if isinstance(reacting_mol, Sidechain) else reacting_mol._id
            reacting_mol = {'_id': _id, 'kekule': reacting_mol.kekule}
        return cls(rxn_type, AllChem.ReactionFromSmarts(smarts).ToBinary(), smarts, template._id, reacting_mol, rxn_atom_idx)

    @classmethod
    def from_dict(cls, data, _id=None):
        return cls(data['type'], data['binary'], data['smarts'], data['template'], data['reacting_mol'], data['rxn_atom_idx'], _id=_id)

    @property
    def rxn(self):
        return AllChem.ChemicalReaction(self.binary)

    def to_dict(self):
        data = deepcopy(self.__dict__)
        data.pop('_id')
        return data

    def to_reduced_dict(self):
        return {'_id': self._id, 'type': self.type}


class AbstractPrediction:
    def __init__(self, predictions, reacting_mol, solvent, _id):
        self._id = _id
        self.predictions = predictions
        self.reacting_mol = reacting_mol
        self.solvent = solvent

    def to_dict(self):
        data = deepcopy(self.__dict__)
        data.pop('_id')
        return data


class RegioSQMPrediction(AbstractPrediction):
    def __init__(self, predictions, reacting_mol, solvent, cutoff, _id=None):
        super().__init__(predictions, reacting_mol, solvent, _id)
        self.cutoff = cutoff

    @classmethod
    def from_dict(cls, data, _id=None):
        return cls(data['predictions'], data['reacting_mol'], data['solvent'], data['cutoff'], _id=_id)


class pKaPrediction(AbstractPrediction):
    def __init__(self, predictions, reacting_mol, solvent, _id=None):
        super().__init__(predictions, reacting_mol, solvent, _id)

    @classmethod
    def from_dict(cls, data, _id=None):
        return cls(data['predictions'], data['reacting_mol'], data['solvent'], _id=_id)


class PeptidePlan:
    def __init__(self, peptide_length):
        self.combinations = set()
        self.peptide_length = peptide_length
        self.valid_lengths = (peptide_length, peptide_length + 1)

    def __iter__(self):
        return iter(self.combinations)

    def __next__(self):
        return next(self.combinations)

    def __len__(self):
        return len(self.combinations)

    @classmethod
    def from_array(cls, array):
        peptide_plan = cls(len(array[0]))
        peptide_plan.combinations = np.array(array)
        return peptide_plan

    def add(self, combination):
        if isinstance(self.combinations, set) and len(combination) in self.valid_lengths:
            self.combinations.add(combination)
        else:
            raise RuntimeError('Cannot add to peptide plan once uniqueness checks have been switched off!')

    def data(self):
        self.combinations = list(self.combinations)
        self.combinations.sort(key=len)
        return self.combinations
