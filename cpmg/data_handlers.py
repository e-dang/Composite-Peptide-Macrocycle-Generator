from collections import defaultdict
from itertools import product, chain

import cpmg.repository as repo
import cpmg.reactions as rxns
from cpmg.ranges import Key, WholeRange
import cpmg.config as config
import cpmg.generators as generators
import cpmg.utils as utils
from cpmg.models import METHANE


class AbstractDataHandler:
    def __init__(self, saver):
        self.saver = saver

    def save(self, data):
        return self.saver.save(data)


class SidechainDataHandler(AbstractDataHandler):
    STRING = generators.SidechainModifier.STRING

    def __init__(self):
        self.connection_repo = repo.create_connection_repository()
        self.sidechain_repo = repo.create_sidechain_repository()
        super().__init__(self.sidechain_repo)

    def load(self, key=Key(WholeRange())):
        sidechains = self.sidechain_repo.load(key)
        connections = list(self.connection_repo.load(key))
        for sidechain in sidechains:
            yield [sidechain, list(filter(lambda x: x.kekule != sidechain.connection, connections))]


class MonomerDataHandler(AbstractDataHandler):
    STRING = generators.MonomerGenerator.STRING

    def __init__(self):
        self.backbone_repo = repo.create_backbone_repository()
        self.sidechain_repo = repo.create_sidechain_repository()
        super().__init__(repo.create_monomer_repository())

    def load(self, key=Key(WholeRange())):
        sidechains = self.sidechain_repo.load(key)
        backbones = list(self.backbone_repo.load(key))
        for sidechain in sidechains:
            yield [sidechain, backbones]


class PeptidePlanDataHandler(AbstractDataHandler):
    STRING = generators.PeptidePlanGenerator.STRING

    def __init__(self):
        self.monomer_repo = repo.create_monomer_repository()
        super().__init__(repo.create_peptide_plan_repository())

    def load(self, key=Key(WholeRange())):
        return [list(self.monomer_repo.load(key))]


class PeptideDataHandler(AbstractDataHandler):
    STRING = generators.PeptideGenerator.STRING

    def __init__(self):
        self.plan_repo = repo.create_peptide_plan_repository()
        self.monomer_repo = repo.create_monomer_repository()
        super().__init__(repo.create_peptide_repository())

    def load(self, key=Key(WholeRange())):
        self._hash_monomers()

        self.hashed_plan = {}
        for _id, indices in self.plan_repo.load(key):
            monomers = [self.index_hash[index] for index in indices]
            self.hashed_plan[indices] = _id
            yield monomers

    def save(self, data):
        completed_combos = []
        for peptide in data:
            indices = tuple(self.id_hash[monomer['_id']] for monomer in peptide.monomers)
            completed_combos.append(self.hashed_plan[indices])

        self.plan_repo.deactivate_records(Key(completed_combos))
        return self.saver.save(data)

    def _hash_monomers(self):
        self.index_hash = {}
        self.id_hash = {}
        for monomer in self.monomer_repo.load():
            self.index_hash[monomer.index] = monomer
            self.id_hash[monomer._id] = monomer.index


class TemplatePeptideDataHandler(AbstractDataHandler):
    STRING = generators.TemplatePeptideGenerator.STRING

    def __init__(self):
        self.template_repo = repo.create_template_repository()
        self.peptide_repo = repo.create_peptide_repository()
        super().__init__(repo.create_template_peptide_repository())

    def load(self, key=Key(WholeRange())):
        peptides = self.peptide_repo.load()
        templates = list(self.template_repo.load(key))
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
    STRING = generators.InterMolecularReactionGenerator.STRING

    def __init__(self):
        self.sidechain_repo = repo.create_sidechain_repository()
        self.monomer_repo = repo.create_monomer_repository()
        self.template_repo = repo.create_template_repository()
        super().__init__(repo.create_reaction_repository())

    def load(self, key=Key(WholeRange())):
        templates = list(self.template_repo.load(key))
        for component in self._get_filtered_sidechains(key) + self._get_filtered_monomers(key):
            yield [component, templates]

    def _get_filtered_sidechains(self, key):
        return list(filter(lambda x: x.connection == METHANE, self.sidechain_repo.load(key)))

    def _get_filtered_monomers(self, key):
        return list(filter(lambda x: x.required and x.imported, self.monomer_repo.load(key)))


class IntraMolecularReactionDataHandler(AbstractDataHandler):
    STRING = generators.IntraMolecularReactionGenerator.STRING

    def __init__(self):
        self.template_repo = repo.create_template_repository()
        super().__init__(repo.create_reaction_repository())

    def load(self, key=Key(WholeRange())):
        return self.template_repo.load(key)


get_all_handler_strings = utils.get_module_strings(__name__)

create_handler_from_string = utils.create_factory_function_closure(__name__, 'handler')
