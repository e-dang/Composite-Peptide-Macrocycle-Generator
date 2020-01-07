from abc import ABC, abstractmethod
from copy import copy
from random import sample

import utils
import project_io


class IPlanner(ABC):

    @abstractmethod
    def create_plan(self):
        pass


class PeptidePublicationPlanner(IPlanner):

    def __init__(self, monomer_io, peptide_length, num_peptides):
        self.monomers = monomer_io.load()
        self.saver = project_io.PeptidePlannerIO(peptide_length)
        self.peptide_length = peptide_length
        self.num_peptides = num_peptides
        self.monomer_combinations = set()

    def create_plan(self):

        self.validate_num_peptides()
        self.create_minimum_list()
        self.create_remaining_list()
        self.saver.save(self.monomer_combinations)

    def create_minimum_list(self):

        for position in range(self.peptide_length):
            for desired_monomer in self.monomers:
                fillers = self.get_fillers(desired_monomer)
                self.monomer_combinations.add(
                    tuple(fillers[0:position] + [desired_monomer['index']] + fillers[position:]))

    def create_remaining_list(self):

        monomers = [copy(self.monomers) for _ in range(self.peptide_length)]
        while len(self.monomer_combinations) < self.num_peptides:
            random_sample = utils.random_order_cartesian_product(*monomers)
            if self.validate_monomers(random_sample):
                self.monomer_combinations.add(tuple(monomer['index'] for monomer in random_sample))

    def get_fillers(self, desired_monomer):

        while True:
            monomers = list(sample(self.monomers, self.peptide_length - 1))
            monomers.append(desired_monomer)
            if self.validate_monomers(monomers):
                return [monomer['index'] for monomer in monomers[:-1]]

    def validate_monomers(self, monomers):

        if len(list(filter(lambda x: x['required'], monomers))) > 0:
            return True

        return False

    def validate_num_peptides(self):

        if self.num_peptides < len(self.monomers) * self.peptide_length:
            raise ValueError('The requested number of peptides needs to be at least as large as the number of monomers'
                             ' times the peptide length')
