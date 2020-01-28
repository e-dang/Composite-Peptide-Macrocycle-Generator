
import random
from copy import deepcopy
from functools import wraps
from itertools import chain, combinations, product

from rdkit import Chem
from rdkit.Chem import AllChem

import macrocycles.proxies as proxies
import macrocycles.utils as utils


def apply_stereochemistry(original_func):
    """
    Decorator that can be applied to a Generator's generate method in order to apply all permutations of stereochemistry
    to unmarked chiral centers. THe stereochemistry at chiral centers that already have assigned stereochemistry is
    preserved.

    Args:
        original_func (func): The original method to wrap with this decorator.

    Returns:
        list[dict]: A list of molecule documents with stereochemistry applied.
    """

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
    """
    Decorator function that can be applied to a Generator's generate method to apply all permutations of methylation
    patterns to a molecule at any aromatic nitrogen with one hydrogen atom. The original unmethylated molecule is
    preserved and also returned in the list of methylated compounds.

    Args:
        original_func (func): The function to be wrapped with this decorator.

    Returns:
        list[dict]: A list of molecule documents with the methylation patterns applied, as well as the original
            molecule.
    """

    candidate_heteroatoms = Chem.MolFromSmarts('[nH1,NH1]')  # only methylate heterocycle amines
    methyl = Chem.MolFromSmarts('[CH4:1]')
    map_nums = (1, 2)

    @wraps(original_func)
    def methlyate_wrapper(*args, **kwargs):

        data = []
        for original_result in original_func(*args, **kwargs):
            data.append(original_result)

            # find all combinations of candidate methylation sites
            matches = Chem.Mol(original_result['binary']).GetSubstructMatches(candidate_heteroatoms)
            for atom_idx in chain.from_iterable(matches):
                mol = Chem.Mol(original_result['binary'])
                mol.GetAtomWithIdx(atom_idx).SetAtomMapNum(2)
                mol = utils.connect_mols(mol, methyl, map_nums=map_nums)

                # format data
                binary = mol.ToBinary()
                Chem.Kekulize(mol)
                doc = deepcopy(original_result)
                doc['_id'] += 'm' + str(atom_idx)
                doc['binary'] = binary
                doc['kekule'] = Chem.MolToSmiles(mol, kekuleSmiles=True)
                doc['modifications'] += 'm'
                data.append(doc)

            # for i in range(1, len(matches) + 1):
            #     for atom_tuple in combinations(matches, r=i):

            #         # apply methylation
            #         macrocycle = Chem.Mol(original_result['binary'])
            #         atom_maps = []
            #         for atom_map, atom_idx in enumerate(chain.from_iterable(atom_tuple), start=2):
            #             macrocycle.GetAtomWithIdx(atom_idx).SetAtomMapNum(atom_map)
            #             atom_maps.append(atom_map)
            #         for atom in macrocycle.GetAtoms():
            #             if atom.GetAtomMapNum() in atom_maps:
            #                 macrocycle = utils.connect_mols(macrocycle, methyl, map_nums=[1, atom.GetAtomMapNum()])

            #         # format data
            #         binary = macrocycle.ToBinary()
            #         Chem.Kekulize(macrocycle)
            #         doc = deepcopy(original_result)
            #         doc['_id'] += 'm' + 'm'.join(map(str, chain.from_iterable(atom_tuple)))
            #         doc['binary'] = binary
            #         doc['kekule'] = Chem.MolToSmiles(macrocycle, kekuleSmiles=True)
            #         doc['modifications'] += 'm' * i
            #         data.append(doc)

        return data

    return methlyate_wrapper


def carboxyl_to_amide(original_func):
    """
    Decorator function that is used to wrap the MacrocycleGenerator's generate method in order to modify the resulting
    macrocycles by replacing the C-terminus carboxyl group with an amide. The original macrocycle with a carboxyl group
    is preserved and returned with the transformed macrocycle.

    Args:
        original_func (func): The MacrocycleGenerator's generate method.

    Returns:
        list[dict]: A list of molecule documents containing the transformed and untransformed macrocycles.
    """

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
    """
    Decorator function that is used to wrap the BiMolecularReactionGenerator's generate method in order to filter out
    reactions that are not predicted to occur by RegioSQM. This filter only applies RegioSQM predictions to Friedel
    Crafts and Pictet Spangler reactions.

    Args:
        original_func (func): The BiMolecularReactionGenerator's generate method.

    Returns:
        list[dict]: A list of valid reaction documents that are predicted to occur by RegioSQM.
    """

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
    """
    Decorator function that wraps the BiMolecularReactionGenerator's generate method to filter out Tsuji Trost reactions
    that are unlikely to occur as predicted by the pKa of the heteroatom involved in the reaction. If the pKa of the
    heteroatom is higher than the cutoff then the reaction is filtered out.

    Args:
        original_func (func): The BiMOlecularReactionGenerator's generate method.

    Returns:
        list[dict]: A list of valid reaction documents that are predicted to occur based on pKa values.
    """

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


def aldehyde_filter(original_func):
    """
    Decorator function used to wrap the MacrocycleGenerator's generate method in order to filter out macrocycles that
    have an aldehyde present as these functional groups are not likely to be present in drug like compounds.

    Args:
        original_func (func): The MacrocycleGenerator's generate method.

    Returns:
        list[dict]: The list of macrocycle documents where the macrocycle doesn't have an aldehyde present.
    """

    aldehyde = Chem.MolFromSmarts('[CH](=O)')

    @wraps(original_func)
    def aldehyde_filter_wrapper(*args, **kwargs):

        data = []
        for original_result in original_func(*args, **kwargs):
            mol = Chem.Mol(original_result['binary'])
            if not mol.HasSubstructMatch(aldehyde):
                data.append(original_result)

        return data

    return aldehyde_filter_wrapper


def molecular_weight_filter(original_func):

    max_molecular_weight = 1000

    @wraps(original_func)
    def molecular_weight_filter_wrapper(*args, **kwargs):

        data = []
        for original_result in original_func(*args, **kwargs):
            mol = Chem.Mol(original_result['binary'])
            if AllChem.CalcExactMolWt(mol) <= max_molecular_weight:
                data.append(original_result)

        return data

    return molecular_weight_filter_wrapper


def rotatable_bond_filter(original_func):

    max_num_rotatable_bonds = 10

    @wraps(original_func)
    def rotatable_bond_filter_wrapper(*args, **kwargs):

        data = []
        for original_result in original_func(*args, **kwargs):
            mol = Chem.Mol(original_result['binary'])
            if AllChem.CalcNumRotatableBonds(mol) <= max_num_rotatable_bonds:
                data.append(original_result)

        return data

    return rotatable_bond_filter_wrapper


def tpsa_filter(original_func):

    max_tpsa = 200

    @wraps(original_func)
    def tpsa_filter_wrapper(*args, **kwargs):

        data = []
        for original_result in original_func(*args, **kwargs):
            mol = Chem.Mol(original_result['binary'])
            if AllChem.CalcTPSA(mol, includeSandP=True) <= max_tpsa:
                data.append(original_result)

        return data

    return tpsa_filter_wrapper


def attach_c_term_cap(original_func):

    c_term_caps = proxies.SideChainCapProxy()
    carboxyl = Chem.MolFromSmarts('[OH1]C(=O)')
    map_nums = (1, 3)

    @wraps(original_func)
    def c_term_cap_wrapper(*args, **kwargs):

        data = []
        for original_result in original_func(*args, **kwargs):
            macrocycle = Chem.Mol(original_result['binary'])
            data.append(original_result)

            for i, match in enumerate(macrocycle.GetSubstructMatches(carboxyl)):
                mol = deepcopy(macrocycle)
                for atom_idx in match:
                    atom = mol.GetAtomWithIdx(atom_idx)
                    if atom.GetSymbol() == 'O' and atom.GetTotalNumHs() == 1:
                        oxygen_idx = atom.GetIdx()
                    elif atom.GetSymbol() == 'C':
                        atom.SetAtomMapNum(1)

                mol = Chem.RWMol(mol)
                mol.RemoveAtom(oxygen_idx)
                rand_cap = random.randint(0, len(c_term_caps) - 1)
                mol = utils.connect_mols(mol, c_term_caps[rand_cap], map_nums=map_nums)

                binary = mol.ToBinary()
                doc = deepcopy(original_result)
                doc['_id'] += str(i) + 'c'
                doc['binary'] = binary
                doc['kekule'] = Chem.MolToSmiles(mol, kekuleSmiles=True)
                doc['modifications'] += 'c'
                data.append(doc)

        return data

    return c_term_cap_wrapper
