
from functools import wraps
from copy import deepcopy
from itertools import product, chain, combinations

from rdkit import Chem

import utils


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


def methylate(original_func):

    @wraps(original_func)
    def methlyate_wrapper(*args, **kwargs):

        return original_func(*args, **kwargs)

    return methlyate_wrapper


def carboxyl_to_amide(original_func):

    carboxyl = Chem.MolFromSmarts('[OH1]C(=O)')

    @wraps(original_func)
    def c_to_a_wrapper(*args, **kwargs):

        data = []
        for original_result in original_func(*args, **kwargs):
            # find all combinations of carboxyls
            macrocycle = Chem.Mol(original_result['binary'])
            matches = macrocycle.GetSubstructMatches(carboxyl)
            for i in range(1, len(matches) + 1):
                for atom_tuple in combinations(matches, r=i):
                    # perform conversion to amide
                    macrocycle = Chem.Mol(original_result['binary'])
                    for atom_idx in chain.from_iterable(atom_tuple):
                        atom = macrocycle.GetAtomWithIdx(atom_idx)
                        if atom.GetSymbol() == 'O' and atom.GetTotalNumHs() == 1:
                            atom.SetAtomicNum(7)  # nitrogen atomic number
                            Chem.SanitizeMol(macrocycle)

                    # format data
                    binary = macrocycle.ToBinary()
                    Chem.Kekulize(macrocycle)
                    doc = deepcopy(original_result)
                    doc['_id'] += 'a' + str(i)
                    doc['binary'] = binary
                    doc['kekule'] = Chem.MolToSmiles(macrocycle, kekuleSmiles=True)
                    doc['modifications'] += 'a' * i
                    data.append(doc)

            data.append(original_result)  # save original result

        return data

    return c_to_a_wrapper


def regiosqm_filter(original_func):

    predictions = utils.get_hashed_regiosqm_predictions()

    @wraps(original_func)
    def regiosqm_filter_wrapper(*args, **kwargs):

        data = []
        for original_result in original_func(*args, **kwargs):
            if original_result['type'] in ('friedel_crafts', 'pictet_spangler'):  # filter is applicable
                sidechain = original_result['sidechain']
                rxn_atom = original_result['rxn_atom_idx']
                if rxn_atom in predictions[sidechain]:  # the reaction is predicted by RegioSQM
                    data.append(original_result)
            else:  # filter is not applicable to this type of reaction
                data.append(original_result)

        return data

    return regiosqm_filter_wrapper
