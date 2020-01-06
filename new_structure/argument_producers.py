from abc import ABC, abstractmethod
from collections import defaultdict
from itertools import product
from random import sample


class IArgumentProducer(ABC):

    @abstractmethod
    def __call__(self, data):
        pass


class NullArgProducer(IArgumentProducer):

    def __call__(self, data):

        return data


class CartesianProductArgProducer(IArgumentProducer):

    def __call__(self, data):
        return product(*data)


class TestPeptideGeneratorArgProducer(IArgumentProducer):

    def __call__(self, data):
        monomers, peptide_length = data

        i = 0
        while i < 100:
            monomer_tup = list(filter(lambda x: x['required'], sample(monomers, peptide_length)))
            if len(monomer_tup) == data.peptide_length:
                i += 1
                yield monomer_tup


class PeptideGeneratorArgProducer(IArgumentProducer):

    def __init__(self, monomer_io):
        self.monomers = sorted(list(monomer_io.load()), key=lambda x: x['index'])

    def __call__(self, data):

        for monomer_idxs in data:
            yield [self.monomers[monomer_idx - 1] for monomer_idx in monomer_idxs]


class MacrocycleGeneratorArgProducer(IArgumentProducer):

    def __call__(self, data):

        # hash all reactions based on sidechain
        reactions = defaultdict(list)
        template_only_reactions = defaultdict(list)
        for reaction in data.reactions:
            try:
                reactions[reaction['sidechain']].append(reaction)
            except KeyError:  # template only reaction
                template_only_reactions[reaction['template']].append(reaction)

        for template_peptide in data.template_peptides:
            args = []
            sidechain_data = set(monomer['sidechain'] for monomer in template_peptide['peptide']['monomers'] if
                                 monomer['sidechain'] is not None)
            for sidechain_id in sidechain_data:
                args.extend(filter(lambda x: x['template'] in ('all', template_peptide['template']),
                                   reactions[sidechain_id]))

            # can apply pictet_spangler reaction
            if template_peptide['template'] != 'temp1':
                try:
                    first_monomer = template_peptide['peptide']['monomers'][0]['sidechain']
                # first monomer has side_chain without an _id (natural amino acid or modified proline)
                except TypeError:
                    first_monomer = None

                # get pictet_spangler reactions involving only the first monomer or any other type of reaction
                args = filter(lambda x: ((x['type'] == 'pictet_spangler' and x['sidechain'] == first_monomer)
                                         or x['type'] in ('friedel_crafts', 'tsuji_trost', 'pyrrolo_indolene')), args)

                # sort args based on if they are pictet_spangler or not
                pictet, other = [], []
                for rxn in args:
                    other.append(rxn) if rxn['type'] != 'pictet_spangler' else pictet.append(rxn)

                first_set = list(product(pictet, other))
                second_set = list(product(template_only_reactions[template_peptide['template']], other))
                if first_set and second_set:
                    for arg in first_set + second_set:
                        yield template_peptide, arg
                elif second_set:
                    for arg in second_set:
                        yield template_peptide, arg
                elif first_set:
                    for arg in first_set:
                        yield template_peptide, arg
                else:
                    for arg in other:
                        yield template_peptide, [arg]
            else:  # friedel_crafts, tsuji_trost, pyrrolo_indolene
                for arg in args:
                    yield template_peptide, [arg]
