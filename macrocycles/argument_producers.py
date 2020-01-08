from abc import ABC, abstractmethod
from collections import defaultdict
from itertools import product, chain
from random import sample

import project_io


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

    def __init__(self):

        self.monomer_io = project_io.get_monomer_io()
        self.loaded = False

    def __call__(self, data):

        if not self.loaded:
            self.monomers = sorted(list(self.monomer_io.load()), key=lambda x: x['index'])
            self.loaded = True

        for monomer_idxs in data:
            monomer_idxs = map(int, monomer_idxs.strip('\n').split(','))
            yield [self.monomers[monomer_idx - 1] for monomer_idx in monomer_idxs]


class MacrocycleGeneratorArgProducer(IArgumentProducer):

    def __call__(self, data):

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
