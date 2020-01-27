import os
from abc import ABC, abstractmethod
from copy import deepcopy
from random import choice, choices

import macrocycles.project_io as project_io
import macrocycles.utils as utils


class IPlanner(ABC):

    @abstractmethod
    def create_plan(self):
        pass


class PeptidePublicationPlanner(IPlanner):

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
        monomer_tuples = utils.random_order_cartesian_product(*monomers)
        while len(self.monomer_combinations) < self.num_peptides:
            for random_sample in monomer_tuples:
                if self.validate_monomers(random_sample):
                    self.monomer_combinations.add(tuple(monomer['index'] for monomer in random_sample))
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

        if self.peptide_length < 5 and 3 > len(list(filter(lambda x: x['required'], monomers))) > 0:
            return True

        if self.peptide_length == 5 and 4 > len(list(filter(lambda x: x['required'], monomers))) > 0:
            return True

        return False

    def c_cap_eligible(self, monomers):

        if self.peptide_length < 5 and 2 > len(list(filter(lambda x: x['required'], monomers))):
            return True

        if self.peptide_length == 5 and 3 > len(list(filter(lambda x: x['required'], monomers))):
            return True

        return False

    def validate_num_peptides(self):

        if self.num_peptides < len(self.monomers) * self.peptide_length:
            raise ValueError('The requested number of peptides needs to be at least as large as the number of monomers'
                             ' times the peptide length')
