
from copy import deepcopy
from functools import wraps
from itertools import chain, combinations, product

from rdkit import Chem

import cpmg.utils as utils
from cpmg.models import CARBOXYL


def apply_stereochemistry(original_func):
    """
    Decorator that can be applied to a Generator's generate method in order to apply all permutations of stereochemistry
    to unmarked chiral centers. The stereochemistry at chiral centers that already have assigned stereochemistry is
    preserved.
    """

    stereochem_types = [Chem.ChiralType.CHI_TETRAHEDRAL_CCW, Chem.ChiralType.CHI_TETRAHEDRAL_CW]

    @wraps(original_func)
    def stereochemistry_wrapper(*args, **kwargs):

        data = []
        for original_mol in original_func(*args, **kwargs):

            mol = original_mol.mol
            filter_func = lambda x: x[1] == '?'  # rdkit marks unassigned stereocenters as '?'
            stereocenters = list(filter(filter_func, Chem.FindMolChiralCenters(mol, includeUnassigned=True)))

            atom_idx_stereo_pairs = [product([atom_idx], stereochem_types) for atom_idx, _ in stereocenters]
            for stereo_assignments in product(*atom_idx_stereo_pairs):
                new_rdkit_mol = original_mol.mol
                for atom_idx, stereo in stereo_assignments:
                    new_rdkit_mol.GetAtomWithIdx(atom_idx).SetChiralTag(stereo)

                # format data
                new_model_mol = deepcopy(original_mol)
                new_model_mol.binary = new_rdkit_mol.ToBinary()
                Chem.Kekulize(new_rdkit_mol)
                new_model_mol.kekule = Chem.MolToSmiles(new_rdkit_mol, kekuleSmiles=True)
                data.append(new_model_mol)

        return data

    return stereochemistry_wrapper


def methylate(original_func):
    """
    Decorator function that can be applied to a Generator's generate method to apply all permutations of methylation
    patterns to a molecule at any aromatic nitrogen with one hydrogen atom. The original unmethylated molecule is
    preserved and also returned in the list of methylated compounds.
    """

    candidate_heteroatoms = Chem.MolFromSmarts('[nH1,NH1]')  # only methylate heterocycle amines
    METHYL_MAP_NUM = 1
    HETERO_ATOM_MAP_NUM = 2
    methyl = Chem.MolFromSmarts(f'[CH4:{METHYL_MAP_NUM}]')
    map_nums = (METHYL_MAP_NUM, HETERO_ATOM_MAP_NUM)

    @wraps(original_func)
    def methlyate_wrapper(*args, **kwargs):

        data = []
        for original_mol in original_func(*args, **kwargs):
            data.append(original_mol)

            # find all combinations of candidate methylation sites
            matches = original_mol.mol.GetSubstructMatches(candidate_heteroatoms)
            for atom_idx in chain.from_iterable(matches):
                mol = original_mol.mol
                mol.GetAtomWithIdx(atom_idx).SetAtomMapNum(HETERO_ATOM_MAP_NUM)
                mol = utils.connect_mols(mol, methyl, map_nums=map_nums)
                new_mol = type(original_mol).add_modification(original_mol, mol, 'm')
                data.append(new_mol)

        return data

    return methlyate_wrapper


def carboxyl_to_amide(original_func):
    """
    Decorator function that is used to wrap the MacrocycleGenerator's generate method in order to modify the resulting
    macrocycles by replacing the C-terminus carboxyl group with an amide. The original macrocycle with a carboxyl group
    is preserved and returned with the transformed macrocycle.
    """

    @wraps(original_func)
    def c_to_a_wrapper(*args, **kwargs):

        data = []
        for original_mol in original_func(*args, **kwargs):
            data.append(original_mol)

            matches = original_mol.mol.GetSubstructMatches(CARBOXYL)
            for i in range(1, len(matches) + 1):
                for atom_tuple in combinations(matches, r=i):

                    # perform conversion to amide
                    macrocycle = original_mol.mol
                    for atom_idx in chain.from_iterable(atom_tuple):
                        atom = macrocycle.GetAtomWithIdx(atom_idx)
                        if atom.GetSymbol() == 'O' and atom.GetTotalNumHs() == 1:
                            atom.SetAtomicNum(7)  # nitrogen atomic number
                            Chem.SanitizeMol(macrocycle)

                    new_mol = type(original_mol).add_modification(original_mol, macrocycle, 'a')
                    data.append(new_mol)

        return data

    return c_to_a_wrapper
