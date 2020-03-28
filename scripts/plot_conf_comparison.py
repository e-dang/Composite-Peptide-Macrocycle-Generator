from rdkit import Chem

from rdkit.Chem import AllChem, Draw

import matplotlib.pyplot as plt
import macrocycles.config as config
import numpy as np


BOLTZMAN_CONSTANT = 0.0019872041


def get_force_fields(mol):
    mol_props = AllChem.MMFFGetMoleculeProperties(mol)
    mol_props.SetMMFFVariant(config.FORCE_FIELD)
    mol_props.SetMMFFDielectricConstant(config.DIELECTRIC)
    return list(map(lambda x: AllChem.MMFFGetMoleculeForceField(
        mol, mol_props, confId=x, ignoreInterfragInteractions=False), range(mol.GetNumConformers())))


def calc_energies(force_fields):
    return np.array(list(map(lambda x: x.CalcEnergy(), force_fields)))


def boltzman(mol, indicies, temp=293.15):  # temp in kelvin
    energies = calc_energies(get_force_fields(mol))
    specific_energies = energies[indicies]  # cis or trans only energies
    energy_min = np.min(specific_energies)

    temp_indices = set()
    for i, energy in enumerate(energies):
        if energy <= energy_min + 5.0:
            temp_indices.add(i)

    desired_indices = list(temp_indices.intersection(indicies))
    desired_energies = energies[desired_indices]
    boltzman_energies = np.exp(-1.0 * desired_energies / (BOLTZMAN_CONSTANT * temp))
    boltzman_sum = np.sum(boltzman_energies)
    return boltzman_energies / boltzman_sum, desired_indices


def boltzman_get_relative(boltzman_distr):
    return boltzman_distr / max(boltzman_distr)


def draw_mol(mol):
    for i in range(mol.GetNumConformers()):
        mol.RemoveConformer(i)
    Draw.MolToImage(mol, size=(750, 750), includeAtomNumbers=True).show()


def calc_distance(mol, atom1, atom2):
    """
    Calculates the distance between two atoms on a given molecule.

    Args:
        mol (RDKit Mol): The molecule containing the two atoms.
        atom1 (int): The index of the first atom on the molecule.
        atom2 (int): The index if the second atom on the molecule.

    Returns:
        int: The distance between the two atoms.
    """

    atom1_position = mol.GetAtomPosition(atom1)
    atom2_position = mol.GetAtomPosition(atom2)
    return atom1_position.Distance(atom2_position)


def compare_distances(our_avg_dists, reported_dists):
    our_error = abs(our_avg_dists - reported_dists)
    return np.average(our_error)


def check_bounds(our_dists, bounds):
    truth_values = []
    actual_values = []
    for distance, (low_bound, upper_bound) in zip(our_dists, bounds):
        truth_values.append(low_bound <= distance <= upper_bound)
        if not truth_values[-1]:
            actual_values.append(abs(distance - upper_bound))
    # print('actual values', actual_values)
    return ~np.array(truth_values)


class Cyclopeptide:
    def __init__(self, atom_pairs, paper_measurements, nmr_measurements, nmr_bounds):
        self.atom_pairs = atom_pairs
        self.paper_measurements = np.array(paper_measurements)
        self.nmr_measurements = np.array(nmr_measurements)
        self.nmr_bounds = nmr_bounds


class Cyclopeptide1CisDeprotonated(Cyclopeptide):
    def __init__(self):
        cis_atom_pairs = []
        cis_atom_pairs.append((43, 70))  # 0 ser HN ser Hb1
        cis_atom_pairs.append((43, 69))  # 1 ser HN ser Hb2
        cis_atom_pairs.append((43, 50))  # 2 ser HN val Ha
        cis_atom_pairs.append((43, 36))  # 3 ser HN pro Ha
        cis_atom_pairs.append((45, 44))  # 4 leu HN ser Ha
        cis_atom_pairs.append((45, 46))  # 5 leu HN leu Ha
        cis_atom_pairs.append((47, 46))  # 6 asp HN leu Ha
        cis_atom_pairs.append((47, 48))  # 7 asp HN asp Ha
        cis_atom_pairs.append((47, 58))  # 8 asp HN asp hb1
        cis_atom_pairs.append((47, 59))  # 9 asp HN asp hb2
        cis_atom_pairs.append((49, 47))  # 10 val HN asp HN
        cis_atom_pairs.append((49, 48))  # 11 val HN asp Ha
        cis_atom_pairs.append((49, 58))  # 12 val HN asp Hb2
        cis_atom_pairs.append((49, 59))  # 13 val HN asp Hb1
        cis_atom_pairs.append((49, 50))  # 14 val HN val Ha
        cis_atom_pairs.append((49, 51))  # 15 val HN val Hb
        cis_atom_pairs.append((50, 51))  # 16 val Ha val Hb
        cis_atom_pairs.append((50, 36))  # 17 val Ha pro Ha

        paper_cis_measurements = []
        paper_cis_measurements.append(2.5)
        paper_cis_measurements.append(2.8)
        paper_cis_measurements.append(3.6)
        paper_cis_measurements.append(3.1)
        paper_cis_measurements.append(2.2)
        paper_cis_measurements.append(2.9)
        paper_cis_measurements.append(2.1)
        paper_cis_measurements.append(2.9)
        paper_cis_measurements.append(2.5)
        paper_cis_measurements.append(2.5)
        paper_cis_measurements.append(2.3)
        paper_cis_measurements.append(2.6)
        paper_cis_measurements.append(2.7)
        paper_cis_measurements.append(2.3)
        paper_cis_measurements.append(2.8)
        paper_cis_measurements.append(2.7)
        paper_cis_measurements.append(2.5)
        paper_cis_measurements.append(2.1)

        nmr_cis_measurements = []
        nmr_cis_measurements.append(2.6)
        nmr_cis_measurements.append(2.6)
        nmr_cis_measurements.append(2.7)
        nmr_cis_measurements.append(3.4)
        nmr_cis_measurements.append(2.1)
        nmr_cis_measurements.append(2.1)
        nmr_cis_measurements.append(2.0)
        nmr_cis_measurements.append(3.4)
        nmr_cis_measurements.append(2.6)
        nmr_cis_measurements.append(2.6)
        nmr_cis_measurements.append(3.5)
        nmr_cis_measurements.append(3.0)
        nmr_cis_measurements.append(3.0)
        nmr_cis_measurements.append(3.0)
        nmr_cis_measurements.append(2.5)
        nmr_cis_measurements.append(2.5)
        nmr_cis_measurements.append(2.3)
        nmr_cis_measurements.append(1.9)

        nmr_cis_bounds = []
        nmr_cis_bounds.append((1.6, 3.8))
        nmr_cis_bounds.append((1.6, 3.8))
        nmr_cis_bounds.append((1.7, 3.7))
        nmr_cis_bounds.append((1.9, 4.2))
        nmr_cis_bounds.append((1.5, 3.1))
        nmr_cis_bounds.append((1.5, 3.1))
        nmr_cis_bounds.append((1.5, 3.0))
        nmr_cis_bounds.append((2.7, 4.2))
        nmr_cis_bounds.append((1.6, 3.9))
        nmr_cis_bounds.append((1.6, 3.9))
        nmr_cis_bounds.append((2.5, 4.5))
        nmr_cis_bounds.append((2.2, 4.1))
        nmr_cis_bounds.append((1.9, 3.7))
        nmr_cis_bounds.append((1.9, 3.7))
        nmr_cis_bounds.append((1.7, 3.2))
        nmr_cis_bounds.append((1.5, 3.5))
        nmr_cis_bounds.append((1.5, 3.3))
        nmr_cis_bounds.append((1.5, 2.6))

        super().__init__(cis_atom_pairs, paper_cis_measurements, nmr_cis_measurements, nmr_cis_bounds)


class Cyclopeptide1TransDeprotonated(Cyclopeptide):
    def __init__(self):
        trans_atom_pairs = []
        trans_atom_pairs.append((43, 44))  # ser HN ser Ha
        trans_atom_pairs.append((45, 44))  # leu HN ser Ha
        trans_atom_pairs.append((47, 46))  # asp HN leu Ha
        trans_atom_pairs.append((47, 48))  # asp HN asp Ha
        trans_atom_pairs.append((47, 58))  # asp HN asp hb1
        trans_atom_pairs.append((47, 59))  # asp HN asp hb2
        trans_atom_pairs.append((49, 48))  # val HN asp Ha
        trans_atom_pairs.append((49, 58))  # val HN asp Hb1
        trans_atom_pairs.append((49, 59))  # val HN asp Hb2
        trans_atom_pairs.append((49, 51))  # val HN val Hb

        paper_trans_measurements = []
        paper_trans_measurements.append(2.8)
        paper_trans_measurements.append(2.1)
        paper_trans_measurements.append(2.2)
        paper_trans_measurements.append(2.9)
        paper_trans_measurements.append(2.4)
        paper_trans_measurements.append(2.5)
        paper_trans_measurements.append(3.4)
        paper_trans_measurements.append(2.8)
        paper_trans_measurements.append(2.4)
        paper_trans_measurements.append(2.5)

        nmr_trans_measurements = []
        nmr_trans_measurements.append(3.0)
        nmr_trans_measurements.append(2.2)
        nmr_trans_measurements.append(2.7)
        nmr_trans_measurements.append(2.5)
        nmr_trans_measurements.append(2.8)
        nmr_trans_measurements.append(2.8)
        nmr_trans_measurements.append(2.5)
        nmr_trans_measurements.append(2.7)
        nmr_trans_measurements.append(2.7)
        nmr_trans_measurements.append(2.6)

        nmr_trans_bounds = []
        nmr_trans_bounds.append((2.0, 4.0))
        nmr_trans_bounds.append((1.5, 3.2))
        nmr_trans_bounds.append((1.7, 3.7))
        nmr_trans_bounds.append((1.5, 3.5))
        nmr_trans_bounds.append((1.8, 4.2))
        nmr_trans_bounds.append((1.8, 4.2))
        nmr_trans_bounds.append((1.5, 3.5))
        nmr_trans_bounds.append((1.7, 3.8))
        nmr_trans_bounds.append((1.7, 3.8))
        nmr_trans_bounds.append((1.6, 3.6))

        super().__init__(trans_atom_pairs, paper_trans_measurements, nmr_trans_measurements, nmr_trans_bounds)


class Cyclopeptide1Cis(Cyclopeptide):
    def __init__(self):
        cis_atom_pairs = []
        cis_atom_pairs.append((36, 39))  # 0 ser HN ser Hb1
        cis_atom_pairs.append((36, 38))  # 1 ser HN ser Hb2
        cis_atom_pairs.append((36, 58))  # 2 ser HN val Ha
        cis_atom_pairs.append((36, 66))  # 3 ser HN pro Ha
        cis_atom_pairs.append((41, 37))  # 4 leu HN ser Ha
        cis_atom_pairs.append((41, 42))  # 5 leu HN leu Ha
        cis_atom_pairs.append((52, 42))  # 6 asp HN leu Ha
        cis_atom_pairs.append((52, 53))  # 7 asp HN asp Ha
        cis_atom_pairs.append((52, 55))  # 8 asp HN asp hb1
        cis_atom_pairs.append((52, 54))  # 9 asp HN asp hb2
        cis_atom_pairs.append((57, 52))  # 10 val HN asp HN
        cis_atom_pairs.append((57, 53))  # 11 val HN asp Ha
        cis_atom_pairs.append((57, 54))  # 12 val HN asp Hb2
        cis_atom_pairs.append((57, 55))  # 13 val HN asp Hb1
        cis_atom_pairs.append((57, 58))  # 14 val HN val Ha
        cis_atom_pairs.append((57, 59))  # 15 val HN val Hb
        cis_atom_pairs.append((58, 59))  # 16 val Ha val Hb
        cis_atom_pairs.append((58, 66))  # 17 val Ha pro Ha

        paper_cis_measurements = []
        paper_cis_measurements.append(2.5)
        paper_cis_measurements.append(2.8)
        paper_cis_measurements.append(3.6)
        paper_cis_measurements.append(3.1)
        paper_cis_measurements.append(2.2)
        paper_cis_measurements.append(2.9)
        paper_cis_measurements.append(2.1)
        paper_cis_measurements.append(2.9)
        paper_cis_measurements.append(2.5)
        paper_cis_measurements.append(2.5)
        paper_cis_measurements.append(2.3)
        paper_cis_measurements.append(2.6)
        paper_cis_measurements.append(2.7)
        paper_cis_measurements.append(2.3)
        paper_cis_measurements.append(2.8)
        paper_cis_measurements.append(2.7)
        paper_cis_measurements.append(2.5)
        paper_cis_measurements.append(2.1)

        nmr_cis_measurements = []
        nmr_cis_measurements.append(2.6)
        nmr_cis_measurements.append(2.6)
        nmr_cis_measurements.append(2.7)
        nmr_cis_measurements.append(3.4)
        nmr_cis_measurements.append(2.1)
        nmr_cis_measurements.append(2.1)
        nmr_cis_measurements.append(2.0)
        nmr_cis_measurements.append(3.4)
        nmr_cis_measurements.append(2.6)
        nmr_cis_measurements.append(2.6)
        nmr_cis_measurements.append(3.5)
        nmr_cis_measurements.append(3.0)
        nmr_cis_measurements.append(3.0)
        nmr_cis_measurements.append(3.0)
        nmr_cis_measurements.append(2.5)
        nmr_cis_measurements.append(2.5)
        nmr_cis_measurements.append(2.3)
        nmr_cis_measurements.append(1.9)

        nmr_cis_bounds = []
        nmr_cis_bounds.append((1.6, 3.8))
        nmr_cis_bounds.append((1.6, 3.8))
        nmr_cis_bounds.append((1.7, 3.7))
        nmr_cis_bounds.append((1.9, 4.2))
        nmr_cis_bounds.append((1.5, 3.1))
        nmr_cis_bounds.append((1.5, 3.1))
        nmr_cis_bounds.append((1.5, 3.0))
        nmr_cis_bounds.append((2.7, 4.2))
        nmr_cis_bounds.append((1.6, 3.9))
        nmr_cis_bounds.append((1.6, 3.9))
        nmr_cis_bounds.append((2.5, 4.5))
        nmr_cis_bounds.append((2.2, 4.1))
        nmr_cis_bounds.append((1.9, 3.7))
        nmr_cis_bounds.append((1.9, 3.7))
        nmr_cis_bounds.append((1.7, 3.2))
        nmr_cis_bounds.append((1.5, 3.5))
        nmr_cis_bounds.append((1.5, 3.3))
        nmr_cis_bounds.append((1.5, 2.6))

        super().__init__(cis_atom_pairs, paper_cis_measurements, nmr_cis_measurements, nmr_cis_bounds)


class Cyclopeptide1Trans(Cyclopeptide):
    def __init__(self):
        trans_atom_pairs = []
        trans_atom_pairs.append((36, 37))  # ser HN ser Ha
        trans_atom_pairs.append((41, 37))  # leu HN ser Ha
        trans_atom_pairs.append((52, 42))  # asp HN leu Ha
        trans_atom_pairs.append((52, 53))  # asp HN asp Ha
        trans_atom_pairs.append((52, 55))  # asp HN asp hb1
        trans_atom_pairs.append((52, 54))  # asp HN asp hb2
        trans_atom_pairs.append((57, 53))  # val HN asp Ha
        trans_atom_pairs.append((57, 55))  # val HN asp Hb1
        trans_atom_pairs.append((57, 54))  # val HN asp Hb2
        trans_atom_pairs.append((57, 59))  # val HN val Hb

        paper_trans_measurements = []
        paper_trans_measurements.append(2.8)
        paper_trans_measurements.append(2.1)
        paper_trans_measurements.append(2.2)
        paper_trans_measurements.append(2.9)
        paper_trans_measurements.append(2.4)
        paper_trans_measurements.append(2.5)
        paper_trans_measurements.append(3.4)
        paper_trans_measurements.append(2.8)
        paper_trans_measurements.append(2.4)
        paper_trans_measurements.append(2.5)

        nmr_trans_measurements = []
        nmr_trans_measurements.append(3.0)
        nmr_trans_measurements.append(2.2)
        nmr_trans_measurements.append(2.7)
        nmr_trans_measurements.append(2.5)
        nmr_trans_measurements.append(2.8)
        nmr_trans_measurements.append(2.8)
        nmr_trans_measurements.append(2.5)
        nmr_trans_measurements.append(2.7)
        nmr_trans_measurements.append(2.7)
        nmr_trans_measurements.append(2.6)

        nmr_trans_bounds = []
        nmr_trans_bounds.append((2.0, 4.0))
        nmr_trans_bounds.append((1.5, 3.2))
        nmr_trans_bounds.append((1.7, 3.7))
        nmr_trans_bounds.append((1.5, 3.5))
        nmr_trans_bounds.append((1.8, 4.2))
        nmr_trans_bounds.append((1.8, 4.2))
        nmr_trans_bounds.append((1.5, 3.5))
        nmr_trans_bounds.append((1.7, 3.8))
        nmr_trans_bounds.append((1.7, 3.8))
        nmr_trans_bounds.append((1.6, 3.6))

        super().__init__(trans_atom_pairs, paper_trans_measurements, nmr_trans_measurements, nmr_trans_bounds)


class Cyclopeptide2(Cyclopeptide):
    def __init__(self):
        atom_pairs = []
        atom_pairs.append((81, 70))  # arg HN arg Ha
        atom_pairs.append((81, 42))  # arg HN gly HN
        atom_pairs.append((81, 62))  # arg HN val Ha
        atom_pairs.append((42, 70))  # gly HN arg Ha
        atom_pairs.append((45, 46))  # asp HN asp Ha
        atom_pairs.append((45, 50))  # asp HN phe HN
        atom_pairs.append((50, 46))  # phe HN asp Ha
        atom_pairs.append((50, 51))  # phe HN phe Ha
        atom_pairs.append((50, 62))  # phe HN val Ha

        paper_measurements = []
        paper_measurements.append(2.30)
        paper_measurements.append(2.07)
        paper_measurements.append(2.08)
        paper_measurements.append(2.27)
        paper_measurements.append(2.88)
        paper_measurements.append(3.92)
        paper_measurements.append(2.20)
        paper_measurements.append(2.88)
        paper_measurements.append(5.00)

        # paper only gives upper bounds so assume lower is 0
        nmr_bounds = []
        nmr_bounds.append((0.0, 2.78))
        nmr_bounds.append((0.0, 3.74))
        nmr_bounds.append((0.0, 2.65))
        nmr_bounds.append((0.0, 2.94))
        nmr_bounds.append((0.0, 2.90))
        nmr_bounds.append((0.0, 3.68))
        nmr_bounds.append((0.0, 2.74))
        nmr_bounds.append((0.0, 3.50))
        nmr_bounds.append((0.0, 5.19))

        super().__init__(atom_pairs, paper_measurements, None, nmr_bounds)


def measure_inter_atomic_dists_cyclopep0():
    # deprotonated smiles O=C([C@H]1N(CCC1)C2=O)N[C@H](C(N[C@@H](C(N[C@H](C(N[C@H]2C(C)C)=O)CC([O-])=O)=O)CC(C)C)=O)CO
    # smiles: O=C(N[C@@H](CO)C(N[C@H](CC(C)C)C(N[C@@H](CC(O)=O)C(N[C@H]1C(C)C)=O)=O)=O)[C@H]2N(C1=O)CCC2
    mol = Chem.MolFromPDBFile('cpep_0.pdb', removeHs=False)
    # draw_mol(mol)
    # exit()

    cis_confs = []
    trans_confs = []
    cis_indices = []
    trans_indices = []
    for i, conf in enumerate(mol.GetConformers()):
        degree = abs(AllChem.GetDihedralDeg(conf, 2, 3, 7, 19))  # deprotonated
        # degree = AllChem.GetDihedralDeg(conf, 22, 31, 30, 29) # protonated
        # print(degree)
        if degree < 0:
            degree = 360 + degree
        if degree <= 90:
            cis_confs.append(conf.GetId())
            cis_indices.append(i)
        else:
            trans_confs.append(conf.GetId())
            trans_indices.append(i)

    # print(cis_indices)
    # print(trans_indices)

    cis_weights, cis_indices = boltzman(mol, cis_indices)
    cis_weights = boltzman_get_relative(cis_weights)
    trans_weights, trans_indices = boltzman(mol, trans_indices)
    trans_weights = boltzman_get_relative(trans_weights)
    # cis_weights = np.ones_like(cis_indices)
    # trans_weights = np.ones_like(trans_indices)

    print(f'cis: {len(cis_indices)}, trans: {len(trans_indices)}, trans/cis ratio: {len(trans_indices)/len(cis_indices)}', )

    cis_ = Cyclopeptide1CisDeprotonated()
    trans_ = Cyclopeptide1TransDeprotonated()

    cis_data = []
    trans_data = []
    for conf in mol.GetConformers():
        data = []
        if conf.GetId() in cis_confs and conf.GetId() in cis_indices:
            for i, (atom1, atom2) in enumerate(cis_.atom_pairs):
                data.append(calc_distance(conf, atom1, atom2))
            cis_data.append(data)
        elif conf.GetId() in trans_confs and conf.GetId() in trans_indices:
            for atom1, atom2 in trans_.atom_pairs:
                data.append(calc_distance(conf, atom1, atom2))
            trans_data.append(data)

    cis_data = np.array(cis_data).reshape((len(cis_indices), -1))
    trans_data = np.array(trans_data).reshape((len(trans_indices), -1))

    cis_data = np.round(np.average(cis_data, axis=0, weights=cis_weights), decimals=1)
    trans_data = np.round(np.average(trans_data, axis=0, weights=trans_weights), decimals=1)
    # cis_data = cis_data[0, :]
    # trans_data = trans_data[0, :]

    avg_cis_error_paper = compare_distances(cis_data, cis_.paper_measurements)
    avg_trans_error_paper = compare_distances(trans_data, trans_.paper_measurements)
    avg_cis_error_nmr = compare_distances(cis_data, cis_.nmr_measurements)
    avg_trans_error_nmr = compare_distances(trans_data, trans_.nmr_measurements)
    avg_cis_paper_nmr_error = compare_distances(cis_.paper_measurements, cis_.nmr_measurements)
    avg_trans_paper_nmr_error = compare_distances(trans_.paper_measurements, trans_.nmr_measurements)
    cis_bounds_check = check_bounds(cis_data, cis_.nmr_bounds)
    trans_bounds_check = check_bounds(trans_data, trans_.nmr_bounds)

    print('cyclopeptide 1')
    print('avg_cis_error_paper', avg_cis_error_paper)
    print('avg_trans_error_paper', avg_trans_error_paper)
    print('avg_cis_error_nmr', avg_cis_error_nmr)
    print('avg_trans_error_nmr', avg_trans_error_nmr)
    print('avg_cis_paper_nmr_error', avg_cis_paper_nmr_error)
    print('avg_trans_paper_nmr_error', avg_trans_paper_nmr_error)
    print('distances out of nmr bounds (cis):', sum(cis_bounds_check))
    print('distances out of nmr bounds (trans):', sum(trans_bounds_check))
    print()

    return cis_data, trans_data


def measure_inter_atomic_dists_cyclopep1():
    # smiles: O=C(NCC(N[C@@H](CC(O)=O)C(N[C@H](CC1=CC=CC=C1)C(N(C)[C@H]2C(C)C)=O)=O)=O)[C@H](CCCNC(N)=N)NC2=O
    mol = Chem.MolFromPDBFile('cyclopeptide_1.pdb', removeHs=False)

    cyclopep = Cyclopeptide2()

    data = []
    for conf in mol.GetConformers():
        conf_data = []
        for atom1, atom2 in cyclopep.atom_pairs:
            conf_data.append(calc_distance(conf, atom1, atom2))
        data.append(conf_data)

    data = np.array(data).reshape((mol.GetNumConformers(), -1))
    avg_distances = np.round(np.average(data, axis=0), decimals=2)
    print(avg_distances)
    avg_error = compare_distances(avg_distances, cyclopep.paper_measurements)
    bounds_check = check_bounds(avg_distances, cyclopep.nmr_bounds)

    print('cyclopeptide2 - cilengitide')
    print('avg_error', avg_error)
    print('distances out of nmr bounds:', sum(bounds_check))

    return avg_distances


def measure_inter_atomic_dists_cyclopep2():
    # smiles: O=C(N[C@@H](CCCNC(N)=N)C(N[C@@H](CC1=CNC2=CC=CC=C12)C(N[C@@H](CC3=CNC4=CC=CC=C34)C(N[C@@H](CCCNC(N)=N)C(N[C@H]5CC6=CC=CC=C6)=O)=O)=O)=O)[C@H](CCCNC(N)=N)NC5=O
    mol = Chem.MolFromPDBFile('cyclopeptide_2.pdb', removeHs=False)


def plot_cyclopep_data(our_avg_distances, cyclopep, name):
    low_bound, upper_bound = zip(*cyclopep.nmr_bounds)
    x = range(1, len(low_bound) + 1)
    fig, ax = plt.subplots(1, 1)
    ax.set_facecolor("silver")
    ax.plot(x, low_bound, c='silver')
    ax.plot(x, upper_bound, c='silver')
    ax.fill_between(x, low_bound, upper_bound, color='white')
    # ax.plot(x, our_avg_distances, label='ConfBuster++', color='mediumblue', marker='+', ms=6.5, lw=0.2, alpha=0.5)
    ax.plot(x, our_avg_distances, color='mediumblue', lw=0.3, alpha=0.5)
    ax.scatter(x, our_avg_distances, label='ConfBuster++', color='mediumblue', marker='+', s=50.0)
    # ax.plot(x, cyclopep.paper_measurements, label='Paper', color='limegreen', marker='o', ms=3.5, lw=0.2, alpha=0.5)
    ax.plot(x, cyclopep.paper_measurements, color='forestgreen', lw=0.3, alpha=0.5)
    ax.scatter(x, cyclopep.paper_measurements, label='Kamenik et al.', color='forestgreen', marker='o', s=25.0)
    ax.plot(x, cyclopep.nmr_measurements, label='NMR', color='red', lw=1.5, alpha=1)
    ax.set_xlim((min(x), max(x)))
    ax.set_ylim((1, max([max(our_avg_distances), max(cyclopep.paper_measurements)]) + 1))
    ax.set_xticks(x)
    ax.set(xlabel='Proton Pair', ylabel='Inter-Proton Distance (Å)')
    ax.legend(loc='upper left')
    plt.savefig(f'comparison_{name}.png', dpi=2000)
    # plt.show()


def plot_cyclopep_data_bar(our_avg_distances, cyclopep, name):
    x = range(1, len(cyclopep.paper_measurements) + 1)
    fig, (ax1, ax2) = plt.subplots(2, 1)
    our_difference = our_avg_distances - cyclopep.nmr_measurements
    paper_difference = cyclopep.paper_measurements - cyclopep.nmr_measurements
    markerline_1, *_ = ax1.stem(x, our_difference, use_line_collection=True,
                                label='ConfBuster++', linefmt='mediumblue', markerfmt='o')
    markerline_2, *_ = ax2.stem(x, paper_difference, use_line_collection=True,
                                label='Kamenik et al.', linefmt='forestgreen', markerfmt='o')
    markerline_1.set_markerfacecolor('mediumblue')
    markerline_1.set_markeredgecolor('mediumblue')
    markerline_1.set_markersize(4)
    markerline_2.set_markerfacecolor('forestgreen')
    markerline_2.set_markeredgecolor('forestgreen')
    markerline_2.set_markersize(4)
    # ax1.set_xlim((min(x) - 0.5, max(x) + 0.5))
    ax2.set_xlim((min(x) - 0.5, max(x) + 0.5))
    ax1.set_ylim((-1.5, 1.5))
    ax2.set_ylim((-1.5, 1.5))
    ax1.set_xticks(x)
    ax2.set_xticks(x)
    fig.legend()
    fig.text(0.04, 0.5, 'Difference From NMR Distances (Å)', ha='center', va='center', rotation='vertical')
    ax2.set(xlabel='Proton Pair')
    plt.savefig(f'differences_{name}.png', dpi=2000)
    # plt.show()


cis_avg_distances, trans_avg_distances = measure_inter_atomic_dists_cyclopep0()  # proline containing compound
avg_distances = measure_inter_atomic_dists_cyclopep1()  # cilengitide
# measure_inter_atomic_dists_cyclopep2()
# exit()
plot_cyclopep_data(cis_avg_distances, Cyclopeptide1Cis(), 'cis')
plot_cyclopep_data(trans_avg_distances, Cyclopeptide1Trans(), 'trans')
# plot_cyclopep_data(avg_distances, Cyclopeptide2(), 'cyclopep2')
# plot_cyclopep_data_bar(cis_avg_distances, Cyclopeptide1Cis(), 'cis')
# plot_cyclopep_data_bar(trans_avg_distances, Cyclopeptide1Trans(), 'trans')
# plot_cyclopep_data_bar(avg_distances, Cyclopeptide2(), 'cyclopep2')
