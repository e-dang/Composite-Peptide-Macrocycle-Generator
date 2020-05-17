from functools import wraps

from rdkit import Chem
from rdkit.Chem import AllChem

import cpmg.config as config
import cpmg.proxies as proxies
import cpmg.reactions as rxns


def regiosqm_filter(original_func):
    predictions = proxies.HashedRegioSQMProxy()
    types = (rxns.FriedelCrafts.TYPE, rxns.PictetSpangler.TYPE)

    @wraps(original_func)
    def regiosqm_filter_wrapper(*args, **kwargs):

        data = []
        for reaction in original_func(*args, **kwargs):
            if reaction.type in types:
                if reaction.rxn_atom_idx in predictions[reaction.reacting_mol['kekule']]:
                    data.append(reaction)
            else:  # filter is not applicable to this type of reaction
                data.append(reaction)

        return data

    return regiosqm_filter_wrapper


def pka_filter(original_func):

    predictions = proxies.HashedpKaProxy()
    cutoff = config.PKA_CUTOFF
    types = (rxns.TsujiTrost.TYPE, )

    @wraps(original_func)
    def pka_filter_wrapper(*args, **kwargs):

        data = []
        for reaction in original_func(*args, **kwargs):
            if reaction.type in types:
                pkas = predictions[reaction.reacting_mol['kekule']]
                rxn_atom = reaction.rxn_atom_idx
                if pkas[str(rxn_atom)] < cutoff and pkas[str(rxn_atom)] > 0:
                    data.append(reaction)
            else:   # filter is not applicable to this type of reaction
                data.append(reaction)

        return data

    return pka_filter_wrapper


def aldehyde_filter(original_func):

    aldehyde = Chem.MolFromSmarts('[CH](=O)')

    @wraps(original_func)
    def aldehyde_filter_wrapper(*args, **kwargs):

        data = []
        for original_mol in original_func(*args, **kwargs):
            if not original_mol.mol.HasSubstructMatch(aldehyde):
                data.append(original_mol)

        return data

    return aldehyde_filter_wrapper


def molecular_weight_filter(original_func):

    @wraps(original_func)
    def molecular_weight_filter_wrapper(*args, **kwargs):

        data = []
        for original_mol in original_func(*args, **kwargs):
            if AllChem.CalcExactMolWt(original_mol.mol) <= config.MAX_MW:
                data.append(original_mol)

        return data

    return molecular_weight_filter_wrapper


def rotatable_bond_filter(original_func):

    @wraps(original_func)
    def rotatable_bond_filter_wrapper(*args, **kwargs):

        data = []
        for original_mol in original_func(*args, **kwargs):
            if AllChem.CalcNumRotatableBonds(original_mol.mol) <= config.MAX_ROTATABLE_BONDS:
                data.append(original_mol)

        return data

    return rotatable_bond_filter_wrapper


def tpsa_filter(original_func):

    @wraps(original_func)
    def tpsa_filter_wrapper(*args, **kwargs):

        data = []
        for original_mol in original_func(*args, **kwargs):
            if AllChem.CalcTPSA(original_mol.mol, includeSandP=True) <= config.MAX_TPSA:
                data.append(original_mol)

        return data

    return tpsa_filter_wrapper
