from collections import defaultdict
from itertools import product, chain

import cpmg.repository as repo
import cpmg.reactions as rxns
import cpmg.ranges as ranges
import cpmg.config as config


class AbstractDataHandler:
    def __init__(self, saver):
        self.saver = saver
        self.chunk = []
        self.saved_ids = []

    def save(self, data):
        return self.saver.save(data)

    # def __del__(self):
    #     self.flush()

    # def save(self, data):
    #     self.chunk.extend(data)
    #     if len(self.chunk) > config.CAPACITY:
    #         self.flush()

    # def flush(self):
    #     if len(self.chunk) > 0:
    #         self.saved_ids.extend(self.saver.save(self.chunk))
    #         print(f'saved {len(self.chunk)} documents')
    #         self.chunk = []


class SidechainDataHandler(AbstractDataHandler):
    def __init__(self):
        self.connection_repo = repo.create_connection_repository()
        self.sidechain_repo = repo.create_sidechain_repository()
        super().__init__(self.sidechain_repo)

    def load(self, key=ranges.WholeRange()):
        sidechains = self.sidechain_repo.load(key)
        connections = list(self.connection_repo.load(key))
        for sidechain in sidechains:
            yield [sidechain, list(filter(lambda x: x.kekule != sidechain.connection, connections))]


class MonomerDataHandler(AbstractDataHandler):
    def __init__(self):
        self.backbone_repo = repo.create_backbone_repository()
        self.sidechain_repo = repo.create_sidechain_repository()
        super().__init__(repo.create_monomer_repository())

    def load(self):
        sidechains = self.sidechain_repo.load()
        backbones = list(self.backbone_repo.load())
        for sidechain in sidechains:
            yield [sidechain, backbones]


# class PeptideDataHandler(AbstractDataHandler):
#     def __init__(self, peptide_length):
#         self.peptide_length = peptide_length
#         self.plan_repo = repo.create_peptide_plan_repository()
#         self.monomer_repo = repo.create_monomer_repository()
#         super().__init__(repo.create_peptide_repository())

#     def load(self):
#         self._hash_monomers()
#         self._hash_plan()
#         self.plan = self.plan_repo.load()
#         for indices in self.plan:
#             monomers = [self.monomers[index] for index in indices]
#             yield monomers

#     def save(self, data):
#         completed_combos = models.PeptidePlan(self.peptide_length)
#         for peptide in data:
#             completed_combos.add(peptide.get_indices())
#             self.saver.save(peptide)

#         self.plan_repo.inactivate_records(completed_combos)

#     def _hash_monomers(self):
#         self.monomers = {}
#         for monomer in self.monomer_repo.load():
#             self.monomers[monomer.index] = monomer


class TemplatePeptideDataHandler(AbstractDataHandler):
    def __init__(self):
        self.template_repo = repo.create_template_repository()
        self.peptide_repo = repo.create_peptide_repository()
        super().__init__(repo.create_template_peptide_repository())

    def load(self):
        peptides = self.peptide_repo.load()
        templates = list(self.template_repo.load())
        for peptide in peptides:
            yield [peptide, templates]


# class MacrocycleDataHandler(AbstractDataHandler):
#     PS = rxns.PictetSpangler.TYPE
#     FC = rxns.FriedelCrafts.TYPE
#     TT = rxns.TsujiTrost.TYPE
#     PI = rxns.Pyrroloindoline.TYPE

#     def __init__(self):
#         self.reaction_repo = repo.create_reaction_repository()
#         self.template_peptide_repo = repo.create_template_peptide_repository()
#         super().__init__(repo.create_macrocycle_repository())

#     def load(self):
#         self._hash_reactions(self.reaction_repo.load())
#         return self._match_reactions(self.template_peptide_repo.load())

#     def _match_reactions(self, template_peptides):

#         for template_peptide in template_peptides:
#             reactions = []
#             sidechain_data = set(monomer.sidechain for monomer in template_peptide.peptide.monomers if
#                                  monomer.sidechain is not None)
#             monomer_data = set(
#                 monomer._id for monomer in template_peptide.peptide.monomers if monomer.sidechain is None)
#             for mol_id in chain(sidechain_data, monomer_data):
#                 reactions.extend(filter(lambda x: x.template == template_peptide.template, self.reactions[mol_id]))

#             first_monomer = template_peptide.peptide.monomers[0].sidechain
#             is_proline = template_peptide.peptide.monomers[0].proline

#             # can apply pictet_spangler reactions or template_only_reactions because
#             #   1. Template 2 and 3 has unmasked aldehyde
#             #   2. First monomer has a sidechain which comes from a set of aromatic heterocycles (either the sidechain
#             #      undergo pictet_spangler or the template can undergo pictet_spangler/aldehyde cyclization)
#             #   3. First monomer doesnt contain any aromatic atoms and is not a proline thus the template will be able
#             #      to undergo the template_only_reactions
#             if template_peptide.pictet_spangler_kekule is not None and \
#                     ((first_monomer is not None) or (first_monomer is None and not is_proline)):

#                 # get pictet_spangler reactions involving only the first monomer or any other type of reaction
#                 reactions = filter(lambda x: ((x.type == self.PS and x.reacting_mol == first_monomer)
#                                               or x.type in (self.FC, self.TT, self.PI)), reactions)

#                 # sort reactions based on if they are pictet_spangler or not
#                 pictet, other = [], []
#                 for rxn in reactions:
#                     other.append(rxn) if rxn.type != self.PS else pictet.append(rxn)

#                 if not pictet or first_monomer is None:
#                     rxn_combinations = product(self.template_only_reactions[template_peptide.template], other)
#                 else:
#                     # if pictet_spangler fails because template is not oligomerized to the peptide at the first monomer
#                     # (as may occur with lysine in the peptide) then template_only_reaction will be applied. Otherwise
#                     # pictet_spangler will occur, blocking the template_only_reaction from taking place
#                     rxn_combinations = product(pictet, self.template_only_reactions[template_peptide.template], other)

#                 yield template_peptide, rxn_combinations

#             else:  # friedel_crafts, tsuji_trost, pyrrolo_indolene
#                 yield template_peptide, [[rxn] for rxn in filter(lambda x: x.type != self.PS, rxns)]

#     def _hash_reactions(self, reactions):
#         self.reactions = defaultdict(list)
#         self.template_only_reactions = defaultdict(list)
#         for reaction in reactions:
#             try:
#                 self.reactions[reaction.reacting_mol['_id']].append(reaction)
#             except KeyError:  # template only reaction
#                 self.template_only_reactions[reaction.template].append(reaction)


class InterMolecularReactionDataHandler(AbstractDataHandler):
    def __init__(self):
        self.sidechain_repo = repo.create_sidechain_repository()
        self.monomer_repo = repo.create_monomer_repository()
        super().__init__(repo.create_reaction_repository())

    def load(self):
        return list(self.sidechain_repo.load()) + list(self.monomer_repo.load())


class IntraMolecularReactionDataHandler(AbstractDataHandler):
    def __init__(self):
        self.template_repo = repo.create_template_repository()
        super().__init__(repo.create_reaction_repository())

    def load(self):
        return self.template_repo.load()
