import os
from abc import ABC, abstractmethod
from copy import deepcopy
from random import choice, choices, sample
import json
import glob

from rdkit import Chem
from rdkit.Chem import AllChem

import macrocycles.config as config
import macrocycles.project_io as project_io
import macrocycles.utils as utils


class IPlanner(ABC):

    @abstractmethod
    def create_plan(self):
        pass


class PeptidePublicationPlanner(IPlanner):

    MAX_MW = config.MAX_MW - 258  # 258 is average of template MW

    def __init__(self, peptide_length, num_peptides):
        self.monomers = project_io.get_filtered_monomer_set()
        if not isinstance(self.monomers, list):
            self.monomers = list(self.monomers)
        self.c_cap_monomers = self.get_c_cap_monomers()
        self.saver = project_io.PeptidePlannerIO(peptide_length)
        self.peptide_length = peptide_length
        self.num_peptides = num_peptides
        self.monomer_combinations = set()

    def create_plan(self):

        if not os.path.exists(utils.attach_file_num(self.saver.FILEPATH, self.peptide_length)):
            self.validate_num_peptides()
            self.create_minimum_list()
            self.create_remaining_list()
            self.saver.save(self.monomer_combinations)

    def create_minimum_list(self):

        for position in range(self.peptide_length):
            for desired_monomer in self.monomers:
                for fillers in self.get_fillers(desired_monomer):
                    self.monomer_combinations.add(
                        tuple(fillers[0:position] + [desired_monomer['index']] + fillers[position:]))

    def create_remaining_list(self):
        monomers = [deepcopy(self.monomers) for _ in range(self.peptide_length)]
        while True:
            for random_sample in utils.random_sample_cartesian_product(*monomers, sample_size=self.num_peptides * 5):
                if len(self.monomer_combinations) > self.num_peptides:
                    break
                if self.validate_monomers(random_sample):
                    # if len(self.monomer_combinations) % 1000 == 0:
                    #     print(len(self.monomer_combinations))
                    self.monomer_combinations.add(tuple(monomer['index'] for monomer in random_sample))
                    if self.c_cap_eligible(random_sample):
                        random_sample += [choice(self.c_cap_monomers)]
                        self.monomer_combinations.add(tuple(monomer['index'] for monomer in random_sample))
            else:
                continue
            break

    def get_fillers(self, desired_monomer):

        while True:
            monomers = list(choices(self.monomers, k=self.peptide_length - 1))
            monomers.append(desired_monomer)
            if self.validate_monomers(monomers):
                yield [monomer['index'] for monomer in monomers[:-1]]
                if self.c_cap_eligible(monomers):
                    c_cap = choice(self.c_cap_monomers)
                    yield [monomer['index'] for monomer in monomers[:-1]] + [c_cap['index']]
                break

    def get_c_cap_monomers(self):

        return list(filter(lambda x: x['backbone'] == 'alpha' and x['connection'] == 'methyl' and x['required'], self.monomers))

    def validate_monomers(self, monomers):

        mw = sum(map(AllChem.CalcExactMolWt, map(lambda x: Chem.Mol(x['binary']), monomers)))
        if mw > PeptidePublicationPlanner.MAX_MW:
            return False

        if self.peptide_length < 5 and 3 > len(list(filter(lambda x: x['required'], monomers))):
            return True

        if self.peptide_length == 5 and 4 > len(list(filter(lambda x: x['required'], monomers))) > 0:
            return True

        return False

    def c_cap_eligible(self, monomers):

        if self.peptide_length < 5 and 2 > len(list(filter(lambda x: x['required'], monomers))):
            return True

        return False

    def validate_num_peptides(self):

        if self.num_peptides < len(self.monomers) * self.peptide_length:
            raise ValueError('The requested number of peptides needs to be at least as large as the number of monomers'
                             ' times the peptide length')


# class ConformerPublicationPlanner(IPlanner):

#     def __init__(self, peptide_length, num_conformers, num_macrocycles):
#         if num_macrocycles is None:
#             self.macrocycle_loader = project_io.get_macrocycle_io(peptide_length=peptide_length, job_num=None)
#         else:
#             self.macrocycle_loader = None
#             self.num_macrocycles = num_macrocycles
#         self.saver = project_io.ConformerPlannerIO(peptide_length)
#         self.peptide_length = peptide_length
#         self.num_conformers = num_conformers

#     def create_plan(self):

#         if not os.path.exists(utils.attach_file_num(self.saver.FILEPATH, self.peptide_length)):
#             count = self.count_macrocycles()
#             macrocycle_idxs = list(sample(range(count), self.num_conformers))
#             macrocycle_idxs.sort()
#             self.saver.save(macrocycle_idxs)

#     def count_macrocycles(self):

#         # total number of macrocycles was given to us
#         if self.macrocycle_loader is None:
#             return self.num_macrocycles

#         for i, _ in enumerate(self.macrocycle_loader.iterate()):
#             pass

#         return i


class ConformerPublicationPlanner(IPlanner):

    def __init__(self, peptide_length, num_conformers, num_macrocycles):
        self.macrocycle_loader = project_io.get_macrocycle_io(peptide_length=peptide_length, job_num=None)
        self.saver = project_io.ConformerPlannerIO(peptide_length)
        self.peptide_length = peptide_length
        self.num_conformers = num_conformers

    def create_plan(self):

        if not os.path.exists(utils.attach_file_num(self.saver.FILEPATH, self.peptide_length)):
            filepaths = glob.glob(utils.attach_file_num(self.macrocycle_loader.FILEPATH, self.peptide_length, '*'))
            count = 0
            selected_filepaths = set()
            while count < self.num_conformers:
                filepath = choice(filepaths)
                if filepath not in selected_filepaths:
                    count += self.count_macrocycles(filepath)
                    selected_filepaths.add(filepath)

            self.saver.save(selected_filepaths)

    def count_macrocycles(self, filepath):
        with open(filepath, 'r') as file:
            return len(list(json.load(file)))
