from abc import ABC, abstractmethod
from random import sample


class IPlanner(ABC):

    @abstractmethod
    def create_plan(self):
        pass


class PeptidePublicationPlanner(IPlanner):

    def __init__(self, monomer_io, peptide_length, num_peptides):
        self.monomers = monomer_io.load()
        self.peptide_length = peptide_length
        self.num_peptides = num_peptides

        if num_peptides < len(self.monomers) * peptide_length:
            raise ValueError('The requested number of peptides needs to be at least as large as the number of monomers'
                             ' times the peptide length')

    def create_plan(self):

        for position in range(self.peptide_length):
            for desired_monomer in self.monomers:
                fillers = self.get_fillers(desired_monomer)
                yield fillers[0:position] + [desired_monomer] + fillers[position:]

    def get_fillers(self, desired_monomer):

        while True:
            fillers = list(sample(self.monomers, self.peptide_length - 1))
            fillers.append(desired_monomer)
            fillers = list(filter(lambda x: x['required'], fillers))
            if len(fillers) == self.peptide_length:
                return fillers[:-1]
