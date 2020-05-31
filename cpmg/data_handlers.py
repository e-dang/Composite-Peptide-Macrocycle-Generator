from collections import defaultdict, namedtuple
from itertools import product, chain

import cpmg.repository as repo
import cpmg.reactions as rxns
from cpmg.ranges import Key, WholeRange
import cpmg.generators as generators
import cpmg.utils as utils
from cpmg.models import METHANE


SidechainDataHandlerTuple = namedtuple('SidechainDataHandlerTuple', 'sidechain connections')
MonomerDataHandlerTuple = namedtuple('MonomerDataHandlerTuple', 'sidechain backbones')
PeptidePlanDataHandlerTuple = namedtuple('PeptidePlanDataHandlerTuple', 'monomers peptide_length num_peptides')
PeptideDataHandlerTuple = namedtuple('PeptideDataHandlerTuple', 'monomers peptide_length')
TemplatePeptideDataHandlerTuple = namedtuple('TemplatePeptideDataHandlerTuple', 'peptide templates')
MacrocycleDataHandlerTuple = namedtuple('MacrocycleDataHandlerTuple', 'template_peptide reaction_combos')
ConformerDataHandlerTuple = namedtuple('ConformerDataHandlerTuple', 'macrocycle')
InterMolecularReactionDataHandlerTuple = namedtuple('InterMolecularReactionDataHandlerTuple', 'nucleophile templates')
IntraMolecularReactionDataHandlerTuple = namedtuple('IntraMolecularReactionDataHandlerTuple', 'reacting_mol')


class AbstractDataHandler:
    def __init__(self, saver):
        self.saver = saver

    def save(self, data):
        return self.saver.save(data)

    def estimate_num_records(self):
        pass


class SidechainDataHandler(AbstractDataHandler):
    STRING = generators.SidechainModifier.STRING

    def __init__(self):
        self.connection_repo = repo.create_connection_repository()
        self.sidechain_repo = repo.create_sidechain_repository()
        super().__init__(self.sidechain_repo)

    def load(self, *, sidechain_key=Key(WholeRange()), connection_key=Key(WholeRange())):
        connections = list(self.connection_repo.load(connection_key))
        for sidechain in self.sidechain_repo.load(sidechain_key):
            yield SidechainDataHandlerTuple(sidechain, list(filter(lambda x: x.kekule != sidechain.connection, connections)))

    def estimate_num_records(self):
        return self.sidechain_repo.get_num_records()


class MonomerDataHandler(AbstractDataHandler):
    STRING = generators.MonomerGenerator.STRING

    def __init__(self):
        self.backbone_repo = repo.create_backbone_repository()
        self.sidechain_repo = repo.create_sidechain_repository()
        super().__init__(repo.create_monomer_repository())

    def load(self, sidechain_key=Key(WholeRange()), backbone_key=Key(WholeRange())):
        backbones = list(self.backbone_repo.load(backbone_key))
        for sidechain in self.sidechain_repo.load(sidechain_key):
            yield MonomerDataHandlerTuple(sidechain, backbones)

    def estimate_num_records(self):
        return self.sidechain_repo.get_num_records()


class PeptidePlanDataHandler(AbstractDataHandler):
    STRING = generators.PeptidePlanGenerator.STRING

    def __init__(self):
        self.monomer_repo = repo.create_monomer_repository()
        super().__init__(repo.create_peptide_plan_repository())

    def load(self, *, peptide_length, num_records, monomer_key=Key(WholeRange())):
        yield PeptidePlanDataHandlerTuple(list(self.monomer_repo.load(monomer_key)), peptide_length, num_records)

    def estimate_num_records(self):
        return self.monomer_repo.get_num_records()


class PeptideDataHandler(AbstractDataHandler):
    STRING = generators.PeptideGenerator.STRING

    def __init__(self):
        self.plan_repo = repo.create_peptide_plan_repository()
        self.monomer_repo = repo.create_monomer_repository()
        super().__init__(repo.create_peptide_repository())

    def load(self, *, peptide_plan_key, monomer_key=Key(WholeRange())):
        self._hash_monomers(monomer_key)

        self.hashed_plan = {}
        for _id, indices in self.plan_repo.load(peptide_plan_key):
            monomers = [self.index_hash[index] for index in indices]
            self.hashed_plan[indices] = _id
            yield PeptideDataHandlerTuple(monomers, peptide_plan_key.peptide_length)

    def save(self, data):
        completed_combos = []
        for peptide in data:
            indices = tuple(self.id_hash[monomer['_id']] for monomer in peptide.monomers)
            completed_combos.append(self.hashed_plan[indices])

        self.plan_repo.mark_complete(Key(completed_combos))
        return self.saver.save(data)

    def estimate_num_records(self):
        return self.plan_repo.get_num_records()

    def _hash_monomers(self, key):
        self.index_hash = {}
        self.id_hash = {}
        for monomer in self.monomer_repo.load(key):
            self.index_hash[monomer.index] = monomer
            self.id_hash[monomer._id] = monomer.index


class TemplatePeptideDataHandler(AbstractDataHandler):
    STRING = generators.TemplatePeptideGenerator.STRING

    def __init__(self):
        self.template_repo = repo.create_template_repository()
        self.peptide_repo = repo.create_peptide_repository()
        super().__init__(repo.create_template_peptide_repository())

    def load(self, *, peptide_key, template_key=Key(WholeRange())):
        templates = list(self.template_repo.load(template_key))
        for peptide in list(self.peptide_repo.load(peptide_key)):
            yield TemplatePeptideDataHandlerTuple(peptide, templates)

    def save(self, data):
        completed_peptides = set(template_peptide.peptide['_id'] for template_peptide in data)
        self.peptide_repo.mark_complete(Key(completed_peptides))
        return self.saver.save(data)

    def estimate_num_records(self):
        return self.peptide_repo.get_num_records()


class MacrocycleDataHandler(AbstractDataHandler):
    STRING = generators.MacrocycleGenerator.STRING
    PS = rxns.PictetSpangler.TYPE
    FC = rxns.FriedelCrafts.TYPE
    TT = rxns.TsujiTrost.TYPE
    PI = rxns.Pyrroloindoline.TYPE

    def __init__(self):
        self.reaction_repo = repo.create_reaction_repository()
        self.template_peptide_repo = repo.create_template_peptide_repository()
        super().__init__(repo.create_macrocycle_repository())

    def load(self, *, template_peptide_key, reaction_key=Key(WholeRange())):
        self._hash_reactions(reaction_key)
        return self._match_reactions(template_peptide_key)

    def save(self, data):
        completed_template_peptides = set(macrocycle.template_peptide for macrocycle in data)
        self.template_peptide_repo.mark_complete(Key(completed_template_peptides))
        return self.saver.save(data)

    def estimate_num_records(self):
        return self.template_peptide_repo.get_num_records()

    def _match_reactions(self, key):
        templates = repo.create_template_repository().load()
        ps_capable_templates = tuple(template._id for template in filter(lambda x: x.pictet_spangler_kekule, templates))

        for template_peptide in list(self.template_peptide_repo.load(key)):
            reactions = []
            sidechain_data = set(monomer['sidechain'] for monomer in template_peptide.monomers if
                                 monomer['sidechain'] is not None)
            monomer_data = set(
                monomer['_id'] for monomer in template_peptide.monomers if monomer['sidechain'] is None)
            for mol_id in chain(sidechain_data, monomer_data):
                reactions.extend(filter(lambda x: x.template == template_peptide.template, self.reactions[mol_id]))

            first_monomer = template_peptide.monomers[0]['sidechain']
            is_proline = template_peptide.monomers[0]['proline']

            # can apply pictet_spangler reactions or template_only_reactions because
            #   1. Template 2 and 3 has unmasked aldehyde
            #   2. First monomer has a sidechain which comes from a set of aromatic heterocycles (either the sidechain
            #      undergo pictet_spangler or the template can undergo pictet_spangler/aldehyde cyclization)
            #   3. First monomer doesnt contain any aromatic atoms and is not a proline thus the template will be able
            #      to undergo the template_only_reactions
            if template_peptide.template in ps_capable_templates and \
                    ((first_monomer is not None) or (first_monomer is None and not is_proline)):

                # get pictet_spangler reactions involving only the first monomer or any other type of reaction
                reactions = filter(lambda x: ((x.type == self.PS and x.reacting_mol == first_monomer)
                                              or x.type in (self.FC, self.TT, self.PI)), reactions)

                # sort reactions based on if they are pictet_spangler or not
                pictet, other = [], []
                for rxn in reactions:
                    other.append(rxn) if rxn.type != self.PS else pictet.append(rxn)

                if not pictet or first_monomer is None:
                    rxn_combinations = product(self.template_only_reactions[template_peptide.template], other)
                else:
                    # if pictet_spangler fails because template is not oligomerized to the peptide at the first monomer
                    # (as may occur with lysine in the peptide) then template_only_reaction will be applied. Otherwise
                    # pictet_spangler will occur, blocking the template_only_reaction from taking place
                    rxn_combinations = product(pictet, self.template_only_reactions[template_peptide.template], other)

                yield MacrocycleDataHandlerTuple(template_peptide, rxn_combinations)

            else:  # friedel_crafts, tsuji_trost, pyrrolo_indolene
                yield MacrocycleDataHandlerTuple(template_peptide,
                                                 [[rxn] for rxn in filter(lambda x: x.type != self.PS, reactions)])

    def _hash_reactions(self, key):
        self.reactions = defaultdict(list)
        self.template_only_reactions = defaultdict(list)
        for reaction in self.reaction_repo.load(key):
            try:
                self.reactions[reaction.reacting_mol['_id']].append(reaction)
            except TypeError:  # template only reaction
                self.template_only_reactions[reaction.template].append(reaction)


class ConformerDataHandler(AbstractDataHandler):
    STRING = generators.ConformerGenerator.STRING

    def __init__(self):
        self.macrocycle_repo = repo.create_macrocycle_repository()
        super().__init__(repo.create_conformer_repository())

    def load(self, macrocycle_key):
        for macrocycle in self.macrocycle_repo.load(macrocycle_key):
            yield ConformerDataHandlerTuple(macrocycle)

    def save(self, data):
        ids = self.saver.save(data)
        self.macrocycle_repo.remove_records(ids)
        return ids

    def estimate_num_records(self):
        return self.macrocycle_repo.get_num_records()


class InterMolecularReactionDataHandler(AbstractDataHandler):
    STRING = generators.InterMolecularReactionGenerator.STRING

    def __init__(self):
        self.sidechain_repo = repo.create_sidechain_repository()
        self.monomer_repo = repo.create_monomer_repository()
        self.template_repo = repo.create_template_repository()
        super().__init__(repo.create_reaction_repository())

    def load(self, sidechain_key=Key(WholeRange()), monomer_key=Key(WholeRange()), template_key=Key(WholeRange())):
        templates = list(self.template_repo.load(template_key))
        for nucleophile in self._get_filtered_sidechains(sidechain_key) + self._get_filtered_monomers(monomer_key):
            yield InterMolecularReactionDataHandlerTuple(nucleophile, templates)

    def estimate_num_records(self):
        return self.sidechain_repo.get_num_records() + self.monomer_repo.get_num_records()

    def _get_filtered_sidechains(self, key):
        return list(filter(lambda x: x.connection == METHANE, self.sidechain_repo.load(key)))

    def _get_filtered_monomers(self, key):
        return list(filter(lambda x: x.required and x.imported, self.monomer_repo.load(key)))


class IntraMolecularReactionDataHandler(AbstractDataHandler):
    STRING = generators.IntraMolecularReactionGenerator.STRING

    def __init__(self):
        self.template_repo = repo.create_template_repository()
        super().__init__(repo.create_reaction_repository())

    def load(self, template_key=Key(WholeRange())):
        for reacting_mol in self.template_repo.load(template_key):
            yield IntraMolecularReactionDataHandlerTuple(reacting_mol)

    def estimate_num_records(self):
        return self.template_repo.get_num_records()


get_all_handler_strings = utils.get_module_strings(__name__)

create_handler_from_string = utils.create_factory_function_closure(__name__, 'handler')
