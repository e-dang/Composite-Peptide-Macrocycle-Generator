
from functools import wraps
from copy import deepcopy
from itertools import product, chain, combinations

from rdkit import Chem

import utils
import proxies


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
                doc['_id'] += 's' + str(i)
                doc['binary'] = new_mol.ToBinary()
                Chem.Kekulize(new_mol)
                doc['kekule'] = Chem.MolToSmiles(new_mol, kekuleSmiles=True)
                data.append(doc)

        return data

    return stereochemistry_wrapper


def methylate(original_func):

    candidate_heteroatoms = Chem.MolFromSmarts('[nH1]')  # only methylate heterocycle amines
    methyl = Chem.MolFromSmarts('[CH4:1]')

    @wraps(original_func)
    def methlyate_wrapper(*args, **kwargs):

        data = []
        for original_result in original_func(*args, **kwargs):
            data.append(original_result)

            # find all combinations of candidate methylation sites
            matches = Chem.Mol(original_result['binary']).GetSubstructMatches(candidate_heteroatoms)
            for i in range(1, len(matches) + 1):
                for atom_tuple in combinations(matches, r=i):

                    # apply methylation
                    macrocycle = Chem.Mol(original_result['binary'])
                    atom_maps = []
                    for atom_map, atom_idx in enumerate(chain.from_iterable(atom_tuple), start=2):
                        macrocycle.GetAtomWithIdx(atom_idx).SetAtomMapNum(atom_map)
                        atom_maps.append(atom_map)
                    for atom in macrocycle.GetAtoms():
                        if atom.GetAtomMapNum() in atom_maps:
                            macrocycle = utils.connect_mols(macrocycle, methyl, map_nums=[1, atom.GetAtomMapNum()])

                    # format data
                    binary = macrocycle.ToBinary()
                    Chem.Kekulize(macrocycle)
                    doc = deepcopy(original_result)
                    doc['_id'] += 'm' + 'm'.join(map(str, chain.from_iterable(atom_tuple)))
                    doc['binary'] = binary
                    doc['kekule'] = Chem.MolToSmiles(macrocycle, kekuleSmiles=True)
                    doc['modifications'] += 'm' * i
                    data.append(doc)

        return data

    return methlyate_wrapper


def carboxyl_to_amide(original_func):

    carboxyl = Chem.MolFromSmarts('[OH1]C(=O)')

    @wraps(original_func)
    def c_to_a_wrapper(*args, **kwargs):

        data = []
        for original_result in original_func(*args, **kwargs):
            data.append(original_result)

            # find all combinations of carboxyls
            matches = Chem.Mol(original_result['binary']).GetSubstructMatches(carboxyl)
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
                    doc['_id'] += 'a' + 'a'.join(map(str, chain.from_iterable(atom_tuple)))
                    doc['binary'] = binary
                    doc['kekule'] = Chem.MolToSmiles(macrocycle, kekuleSmiles=True)
                    doc['modifications'] += 'a' * i
                    data.append(doc)

        return data

    return c_to_a_wrapper


def regiosqm_filter(original_func):

    predictions = proxies.RegioSQMProxy()

    @wraps(original_func)
    def regiosqm_filter_wrapper(*args, **kwargs):

        data = []
        for original_result in original_func(*args, **kwargs):
            if original_result['type'] in ('friedel_crafts', 'pictet_spangler'):  # filter is applicable
                reacting_mol = original_result['reacting_mol']
                rxn_atom = original_result['rxn_atom_idx']
                if rxn_atom in predictions[reacting_mol]:  # the reaction is predicted by RegioSQM
                    data.append(original_result)
            else:  # filter is not applicable to this type of reaction
                data.append(original_result)

        return data

    return regiosqm_filter_wrapper


def pka_filter(original_func):

    predictions = proxies.pKaProxy()
    cutoff = 13.5

    @wraps(original_func)
    def pka_filter_wrapper(*args, **kwargs):

        data = []
        for original_result in original_func(*args, **kwargs):
            if original_result['type'] == 'tsuji_trost':  # filter is applicable
                pkas = predictions[original_result['reacting_mol']]
                rxn_atom = original_result['rxn_atom_idx']
                if pkas[str(rxn_atom)] < cutoff and pkas[str(rxn_atom)] > 0:  # reaction is predicted by pKa filter
                    data.append(original_result)
            else:   # filter is not applicable to this type of reaction
                data.append(original_result)

        return data

    return pka_filter_wrapper
