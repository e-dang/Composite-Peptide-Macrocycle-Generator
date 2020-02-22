from abc import ABC, abstractmethod
from collections import defaultdict
from itertools import chain, product

import macrocycles.project_io as project_io


class IArgumentProducer(ABC):
    """
    Interface for classes that define a way to take data loaded from a DataHandler and format it such that it can be
    used by a specific Generator class.
    """

    @abstractmethod
    def __call__(self, data):
        """
        Abstract method used to format the data into the form expected by the Generator class.

        Args:
            data (iterable): The data that has been loaded by a DataHandler.
        """


class NullArgProducer(IArgumentProducer):
    """
    Implementation of IArgumentProducer that does no formatting to the passed in data. This class is used when the data
    loaded by the DataHandler is already in the correct format.
    """

    def __call__(self, data):
        """
        Method that returns the passed in data without operating on it.

        Args:
            data (iterable): The data to be returned.

        Returns:
            iterable: The data.
        """

        return data


class CartesianProductArgProducer(IArgumentProducer):
    """
    Implementation of IArgumentProducer that returns the cartesian product of the passed in data. This class is used
    when the DataHandler returns multiple types of documents that need to be paired together and then passed to the
    Generator.
    """

    def __call__(self, data):
        """
        Method that creates the cartesian product of the passed in data. The data should be an iterable of iterables in
        this case.

        Args:
            data (iterable[iterable]): Multiple types of documents that need to be paired together and passed to the
                Generator.

        Returns:
            iterable: The paired data documents.
        """
        return product(*data)


class PeptideGeneratorArgProducer(IArgumentProducer):
    """
    Implementation of IArgumentProducer that makes tuples of monomers that are to be formed into a peptide, where the
    monomers to be placed in the tuple are defined by the monomer indices written out by the PeptidePlanner class.
    """

    def __init__(self):
        """
        Initializes instance variables self.monomer_io and self.loaded.
        """

        self.monomer_io = project_io.get_monomer_io()
        self.loaded = False

    def __call__(self, data):
        """
        Method that creates tuples of monomer documents according the the passed in monomer indices.

        Args:
            data (iterable[str]): An iterable of strings that contain the monomer indices for a peptide separated by
                commas.

        Yields:
            iterable(list[dict]): An iterable of lists containing monomer documents to be formed into peptides by the
                PeptideGenerator.
        """

        if not self.loaded:
            self.monomers = sorted(list(self.monomer_io.load()), key=lambda x: x['index'])
            self.loaded = True

        for monomer_idxs in data:
            monomer_idxs = map(int, monomer_idxs.strip('\n').split(','))
            yield [self.monomers[monomer_idx - 1] for monomer_idx in monomer_idxs]


class MacrocycleGeneratorArgProducer(IArgumentProducer):
    """
    Implementation of IArgumentProducer that combines template_peptides and reactions together such that the reactions
    coupled with the template_peptides are those that involve sidechains/monomers in the template_peptide.
    """

    def __call__(self, data):
        """
        Method that combines template_peptides with only the applicable reactions for that template_peptide based on
        its monomer and sidechain composition.

        Args:
            data (tuple[iterable, iterable]): A tuple containing the iterables of template_peptides and reactions.

        Yields:
            tuple[dict, iterable[dict]]: A tuple containing the template_peptide document, along with an iterable of
                applicable reaction documents.
        """

        # hash all reactions based on reacting_mol
        reactions = defaultdict(list)
        template_only_reactions = defaultdict(list)
        for reaction in data.reactions:
            try:
                reactions[reaction['reacting_mol']].append(reaction)
            except KeyError:  # template only reaction
                template_only_reactions[reaction['template']].append(reaction)

        for template_peptide in data.template_peptides:
            rxns = []
            sidechain_data = set(monomer['sidechain'] for monomer in template_peptide['peptide']['monomers'] if
                                 monomer['sidechain'] is not None)
            monomer_data = set(monomer['_id'] for monomer in template_peptide['peptide']
                               ['monomers'] if monomer['sidechain'] is None)
            for mol_id in chain(sidechain_data, monomer_data):
                rxns.extend(filter(lambda x: x['template'] in ('all', template_peptide['template']),
                                   reactions[mol_id]))

            first_monomer = template_peptide['peptide']['monomers'][0]['sidechain']
            is_proline = template_peptide['peptide']['monomers'][0]['is_proline']

            # can apply pictet_spangler reactions or template_only_reactions because
            #   1. Template 2 and 3 has unmasked aldehyde
            #   2. First monomer has a sidechain which comes from a set of aromatic heterocycles (either the sidechain
            #      undergo pictet_spangler or the template can undergo pictet_spangler/aldehyde cyclization)
            #   3. First monomer doesnt contain any aromatic atoms and is not a proline thus the template will be able
            #      to undergo the template_only_reactions
            if template_peptide['template'] != 'temp1' and \
                    ((first_monomer is not None) or (first_monomer is None and not is_proline)):

                # get pictet_spangler reactions involving only the first monomer or any other type of reaction
                rxns = filter(lambda x: ((x['type'] == 'pictet_spangler' and x['reacting_mol'] == first_monomer)
                                         or x['type'] in ('friedel_crafts', 'tsuji_trost', 'pyrrolo_indolene')), rxns)

                # sort rxns based on if they are pictet_spangler or not
                pictet, other = [], []
                for rxn in rxns:
                    other.append(rxn) if rxn['type'] != 'pictet_spangler' else pictet.append(rxn)

                if not pictet or first_monomer is None:
                    rxn_combinations = product(template_only_reactions[template_peptide['template']], other)
                else:
                    # if pictet_spangler fails because template is not oligomerized to the peptide at the first monomer
                    # (as may occur with lysine in the peptide) then template_only_reaction will be applied. Otherwise
                    # pictet_spangler will occur, blocking the template_only_reaction from taking place
                    rxn_combinations = product(pictet, template_only_reactions[template_peptide['template']], other)

                yield template_peptide, rxn_combinations

            else:  # friedel_crafts, tsuji_trost, pyrrolo_indolene
                yield template_peptide, [[rxn] for rxn in filter(lambda x: x['type'] != 'pictet_spangler', rxns)]


class ConformerGeneratorArgProducer(IArgumentProducer):
    """
    Implementation of IArgumentProducer that makes tuples of monomers that are to be formed into a peptide, where the
    monomers to be placed in the tuple are defined by the monomer indices written out by the PeptidePlanner class.
    """

    def __init__(self, **kwargs):
        """
        Initializes instance variables self.monomer_io and self.loaded.
        """

        self.macrocycle_io = project_io.get_macrocycle_io(**kwargs)

    def __call__(self, data):
        """
        Method that creates tuples of monomer documents according the the passed in monomer indices.

        Args:
            data (iterable[str]): An iterable of strings that contain the monomer indices for a peptide separated by
                commas.

        Yields:
            iterable(list[dict]): An iterable of lists containing monomer documents to be formed into peptides by the
                PeptideGenerator.
        """

        try:
            macrocycle_idx = self.get_next_idx(data)
            for i, macrocycle in enumerate(self.macrocycle_io):
                if i == macrocycle_idx:
                    yield macrocycle
                    macrocycle_idx = self.get_next_idx(data)
        except IndexError:
            pass

    def get_next_idx(self, data):

        return int(data.popleft().strip('\n'))
