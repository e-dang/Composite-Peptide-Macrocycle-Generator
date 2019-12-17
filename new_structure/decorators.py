
from functools import wraps
from copy import deepcopy
from itertools import product

from rdkit import Chem

# class RegioSQMFilter:

#     def __init__(self, original_func):
#         self.original_func = original_func

#     @wraps(original_func)
#     def __call__(self, *args, **kwargs):
#         return self.original_func(*args, **kwargs)


def apply_stereochemistry(original_func):

    stereochem_types = [Chem.ChiralType.CHI_TETRAHEDRAL_CCW, Chem.ChiralType.CHI_TETRAHEDRAL_CW]

    @wraps(original_func)
    def stereochemistry_wrapper(*args, **kwargs):

        data = []
        for original_result in original_func(*args, **kwargs):

            # find stereocenters
            parent_mol = Chem.Mol(original_result['binary'])
            filter_func = lambda x: x[1] == '?'
            stereocenters = list(filter(filter_func, Chem.FindMolChiralCenters(parent_mol, includeUnassigned=True)))

            # apply stereochemistry
            atom_idx_stereo_pairs = [product([atom_idx], stereochem_types) for atom_idx, _ in stereocenters]
            for i, stereo_assignments in enumerate(product(*atom_idx_stereo_pairs)):
                new_mol = Chem.Mol(original_result['binary'])
                for atom_idx, stereo in stereo_assignments:
                    new_mol.GetAtomWithIdx(atom_idx).SetChiralTag(stereo)

                # format data
                doc = deepcopy(original_result)
                doc['_id'] = doc['_id'] + str(i)
                doc['binary'] = new_mol.ToBinary()
                Chem.Kekulize(new_mol)
                doc['kekule'] = Chem.MolToSmiles(new_mol, kekuleSmiles=True)
                data.append(doc)

        return data

    return stereochemistry_wrapper


def regiosqm_filter(original_func):

    @wraps(original_func)
    def regiosqm_filter_wrapper(*args, **kwargs):
        return original_func(*args, **kwargs)

    return regiosqm_filter_wrapper
