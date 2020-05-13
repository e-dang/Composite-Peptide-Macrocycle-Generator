import os
from copy import deepcopy
from itertools import chain
from random import choice, choices

from confbusterplusplus.factory import ConfBusterFactory
from rdkit import Chem
from rdkit.Chem import AllChem

import cpmg.config as config
import cpmg.decorators as decorators
import cpmg.filters as filters
import cpmg.models as models
import cpmg.reactions as rxns
import cpmg.repository as repo
import cpmg.utils as temp_utils
import macrocycles.utils as utils
from cpmg.exceptions import InvalidMolecule


class SidechainModifier:
    STRING = models.Sidechain.STRING

    def generate(self, parent_sidechain, connections):
        """
        Method for creating new sidechain molecules from other sidechains by changing the methyl connection to a
        different connection type.
        """

        new_sidechains = []

        # replace the designated attachment position with each type of connection
        for connection in connections:
            for sidechain in Chem.ReplaceSubstructs(parent_sidechain.mol, models.SC_ATTACHMENT_POINT, connection.mol):
                new_sidechains.append(models.Sidechain.from_mol(sidechain, connection, parent_sidechain.shared_id))

        return new_sidechains


class MonomerGenerator:
    STRING = models.Monomer.STRING
    MAP_NUMS = (models.Backbone.MAP_NUM, models.Sidechain.MAP_NUM)

    def generate(self, sidechain, backbones):
        """
        Method that takes a sidechain and backbones creates new monomers.
        """

        sidechain_mol = sidechain.mapped_mol

        temp_utils.clear_isotopes(sidechain_mol)

        monomers = []
        for backbone in backbones:
            monomer = utils.connect_mols(sidechain_mol, backbone.mol, map_nums=self.MAP_NUMS)
            monomers.append(models.Monomer.from_mol(monomer, backbone, sidechain))

        return monomers


class PeptidePlanGenerator:
    STRING = models.PeptidePlan.STRING

    ALPHA_BACKBONE = Chem.MolToSmiles(models.ALPHA_BACKBONE)
    MAX_MW = config.MAX_MW - 258  # 258 is average of template MW
    OVER_SAMPLE_FACTOR = 5

    def generate(self, monomers, peptide_length, num_peptides):
        c_cap_monomers = self._get_c_cap_monomers(monomers)
        peptide_plan = models.PeptidePlan(peptide_length)

        self._create_minimum_list(monomers, c_cap_monomers, peptide_plan)
        self._create_remaining_list(monomers, c_cap_monomers, peptide_plan, num_peptides)

        return [peptide_plan]

    def _create_minimum_list(self, monomers, c_cap_monomers, peptide_plan):

        for position in range(peptide_plan.reg_length):
            for selected_monomer in monomers:
                for fillers in self._get_fillers(selected_monomer, monomers, c_cap_monomers, peptide_plan.reg_length):
                    peptide_plan.add(
                        tuple(fillers[0:position] + [selected_monomer.index] + fillers[position:]))

    def _create_remaining_list(self, monomers, c_cap_monomers, peptide_plan, num_peptides):
        monomers = [deepcopy(monomers) for _ in range(peptide_plan.reg_length)]
        while True:
            for random_sample in utils.random_sample_cartesian_product(*monomers, sample_size=num_peptides * self.OVER_SAMPLE_FACTOR):
                if len(peptide_plan) > num_peptides:
                    break
                if self._validate_monomers(random_sample, peptide_plan.reg_length):
                    peptide_plan.add(tuple(monomer.index for monomer in random_sample))
                    if self._is_c_cap_eligible(random_sample, peptide_plan.reg_length):
                        random_sample += [choice(c_cap_monomers)]
                        peptide_plan.add(tuple(monomer.index for monomer in random_sample))
            else:
                continue
            break

    def _get_fillers(self, desired_monomer, monomer_pool, c_cap_monomers, peptide_length):

        while True:
            monomers = list(choices(monomer_pool, k=peptide_length - 1))
            monomers.append(desired_monomer)
            if self._validate_monomers(monomers, peptide_length):
                yield [monomer.index for monomer in monomers[:-1]]
                if self._is_c_cap_eligible(monomers, peptide_length):
                    c_cap = choice(c_cap_monomers)
                    yield [monomer.index for monomer in monomers[:-1]] + [c_cap.index]
                break

    def _get_c_cap_monomers(self, monomers):
        return list(filter(lambda x: x.backbone['kekule'] == self.ALPHA_BACKBONE and x.connection == models.METHANE and x.required, monomers))

    def _validate_monomers(self, monomers, peptide_length):

        mw = sum(map(AllChem.CalcExactMolWt, map(lambda x: x.mol, monomers)))
        if mw > self.MAX_MW:
            return False

        if peptide_length < 5 and 3 > len(list(filter(lambda x: x.required, monomers))):
            return True

        if peptide_length == 5 and 4 > len(list(filter(lambda x: x.required, monomers))) > 0:
            return True

        return False

    def _is_c_cap_eligible(self, monomers, peptide_length):

        if peptide_length < 5 and 2 > len(list(filter(lambda x: x.required, monomers))):
            return True

        return False


class PeptideGenerator:
    STRING = models.Peptide.STRING

    MONOMER_NITROGEN_MAP_NUM = 1
    PEPTIDE_CARBON_MAP_NUM = 2
    MAP_NUMS = (MONOMER_NITROGEN_MAP_NUM, PEPTIDE_CARBON_MAP_NUM)

    # this reaction only works for alpha amino acids by design
    DECARBOXYLATE = AllChem.ReactionFromSmarts('[*:1]NC([*:2])C(=O)[OH]>>[*:1]NC([*:2])')

    def generate(self, monomers, peptide_length):
        """
        Creates a peptide from a list of monomers.
        """

        for i, monomer in enumerate(monomers):

            if i == 0:
                self.initialize_peptide(monomer)
                continue

            monomer_mol = monomer.mol
            backbone_mol = monomer.backbone_mol

            self.tag_monomer_n_term(monomer_mol, backbone_mol)
            self.tag_peptide_c_term()
            self.add_monomer(monomer_mol, backbone_mol)

        has_c_cap = False
        if len(monomers) != peptide_length:  # c-term cap is present
            self.peptide = self.decarboxylate_c_term()
            has_c_cap = True

        return [models.Peptide.from_mol(self.peptide, peptide_length, has_c_cap, monomers)]

    def initialize_peptide(self, monomer):
        self.pep_size = 1
        self.peptide = monomer.mol
        self.backbone_prev = monomer.backbone_mol

    def tag_monomer_n_term(self, monomer, backbone):
        for atom_idx in monomer.GetSubstructMatch(backbone):
            atom = monomer.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() == 'N' and atom.GetTotalNumHs() != 0:
                atom.SetAtomMapNum(self.MONOMER_NITROGEN_MAP_NUM)

    def tag_peptide_c_term(self):

        flag = False
        for pair in self.peptide.GetSubstructMatches(self.backbone_prev):
            for atom_idx in pair:
                atom = self.peptide.GetAtomWithIdx(atom_idx)
                if atom.GetSymbol() == 'O' and atom.GetTotalNumHs() == 1:  # finds protonated oxygen on carboxyl group
                    self.carboxyl_atom = atom_idx
                elif atom.GetSymbol() == 'C' and atom.GetTotalNumHs() == 0 and \
                        atom.GetHybridization() == Chem.HybridizationType.SP2:  # finds carboxyl carbon
                    attachment_atom = atom

                # this conditional statements determines if the carboxyl oxygen and carbon found above are valid
                # since the nitrogen on a peptide chain should only have 0 or 1 hydrogens (0 if monomer is a proline)
                # or have 2 hydrogens if it is the first monomer in the peptide
                elif (self.pep_size > 1 and atom.GetSymbol() == 'N' and atom.GetTotalNumHs() <= 1) \
                        or self.pep_size == 1:
                    flag = True

            if flag:
                attachment_atom.SetAtomMapNum(self.PEPTIDE_CARBON_MAP_NUM)
                break

    def add_monomer(self, monomer, backbone):
        self.peptide = temp_utils.remove_atom(self.peptide, self.carboxyl_atom)
        self.peptide = utils.connect_mols(self.peptide, monomer, map_nums=self.MAP_NUMS)
        self.pep_size += 1
        self.backbone_prev = backbone

    def decarboxylate_c_term(self):
        products = set()
        for product in chain.from_iterable(self.DECARBOXYLATE.RunReactants((self.peptide,))):
            Chem.SanitizeMol(product)
            products.add(Chem.MolToSmiles(product))

        if len(products) != 1:
            raise RuntimeError('Number of peptides after decarobxylation did not result in a single product!')

        return Chem.MolFromSmiles(products.pop())


class TemplatePeptideGenerator:
    STRING = models.TemplatePeptide.STRING

    # any primary amine or proline n-terminus, but no guanidine
    ELIGIBLE_NITROGENS = Chem.MolFromSmarts('[$([NH2]),$([NH;R]);!$([NH2]C(=O)*);!$([NH2]C(=[NH])[NH]*)]')
    PEPTIDE_NITROGEN_MAP_NUM = 2
    MAP_NUMS = (models.Template.OLIGOMERIZATION_MAP_NUM, PEPTIDE_NITROGEN_MAP_NUM)

    def generate(self, peptide, templates):
        """
        Method that takes a peptide molecule and combines it with each type of template molecule to form
        template-peptide oligomers.
        """

        template_peptides = []
        for template in templates:
            peptide_mol = peptide.mol
            for atom_idx in chain.from_iterable(peptide_mol.GetSubstructMatches(self.ELIGIBLE_NITROGENS)):
                atom = peptide_mol.GetAtomWithIdx(atom_idx)

                atom.SetAtomMapNum(self.PEPTIDE_NITROGEN_MAP_NUM)
                template_peptide = utils.connect_mols(template.oligomerization_mol, peptide_mol, map_nums=self.MAP_NUMS)
                atom.SetAtomMapNum(0)

                template_peptides.append(models.TemplatePeptide.from_mol(template_peptide, template, peptide))

        return template_peptides


class MacrocycleGenerator:
    STRING = models.Macrocycle.STRING

    PS = rxns.PictetSpangler.TYPE
    TPS = rxns.TemplatePictetSpangler.TYPE
    ALDH = rxns.AldehydeCyclization.TYPE
    MAX_ATOM_DIFFERENCE = 5

    @decorators.apply_stereochemistry
    @filters.tpsa_filter
    @filters.rotatable_bond_filter
    @filters.molecular_weight_filter
    @decorators.methylate
    @decorators.carboxyl_to_amide
    @filters.aldehyde_filter
    def generate(self, template_peptide, reaction_combos):

        final_macrocycles = []
        num_atoms = len(template_peptide.mol.GetAtoms())
        for reaction_combo in reaction_combos:
            reactants = [template_peptide.mol]
            reactions = [(reaction.rxn, reaction.type) for reaction in reaction_combo]

            successful_rxn = False
            successful_rxns = []
            for i, (rxn, rxn_type) in enumerate(reactions):

                # pictet_spangler worked, don't use template_only_reactions
                if rxn_type in (self.TPS, self.ALDH) and successful_rxn:
                    continue

                macrocycles = {}
                successful_rxn = False
                for reactant in reactants:
                    for macrocycle in chain.from_iterable(rxn.RunReactants((reactant,))):
                        # santization is expected to fail some what often, mostly due to 2 sidechains where one contains
                        # the other as a substructure which causes the failure
                        if Chem.SanitizeMol(macrocycle, catchErrors=True):
                            continue

                        # protect atoms that participated in this reaction from reacting again in subsequent reactions
                        self._set_protected_atoms(macrocycle, rxn)

                        # the number of atoms shouldn't have changed significantly
                        if abs(num_atoms - len(macrocycle.GetAtoms())) < self.MAX_ATOM_DIFFERENCE:
                            macrocycles[Chem.MolToSmiles(macrocycle)] = macrocycle

                if rxn_type in (self.PS, self.TPS, self.ALDH) and len(macrocycles) == 0:  # reaction failed; reuse reactant
                    continue

                reactants = macrocycles.values()
                successful_rxn = True
                successful_rxns.append(reaction_combo[i])

            for macrocycle in macrocycles:
                try:
                    final_macrocycles.append(models.Macrocycle.from_mol(
                        Chem.MolFromSmiles(macrocycle), '', template_peptide, successful_rxns))
                except InvalidMolecule:
                    continue

        return final_macrocycles

    def _set_protected_atoms(self, macrocycle, reaction):
        for atom in chain.from_iterable(reaction.GetReactingAtoms()):
            atom = macrocycle.GetAtomWithIdx(atom)
            if atom.GetIsAromatic():
                atom.SetProp('_protected', '1')


class ConformerGenerator:
    STRING = models.Conformer.STRING
    ARGS = config.CONFORMER_ARGS

    def __init__(self):
        self.factory = ConfBusterFactory(**self.ARGS._asdict())
        self.factory.MOL_FILE = os.path.join(config.TMP_DIR, 'conf_macrocycle.sdf')
        self.factory.GENETIC_FILE = os.path.join(config.TMP_DIR, 'genetic_results.sdf')
        self.generator = self.factory.create_conformer_generator()

    def generate(self, macrocycle):

        conformer, energies, rmsd, ring_rmsd, *_ = self.generator.generate(Chem.MolFromSmiles(macrocycle.kekule))

        return [models.Conformer.from_macrocycle(conformer, macrocycle, energies, rmsd, ring_rmsd)]


class InterMolecularReactionGenerator:
    STRING = 'inter_' + models.Reaction.STRING

    def __init__(self, impl=None):
        self.impl = impl or rxns.create_intermolecular_reactions()
        self.impl = temp_utils.to_list(self.impl)
        self.backbones = list(map(lambda x: x.mol, repo.create_backbone_repository().load()))

    @filters.pka_filter
    @filters.regiosqm_filter
    def generate(self, nucleophile, templates):
        nucleophile_mol = nucleophile.mol

        non_symmetric_atom_idxs = self._get_non_symmetric_atoms(nucleophile_mol)

        reactions = []
        for atom in [nucleophile_mol.GetAtomWithIdx(atom_idx) for atom_idx in non_symmetric_atom_idxs]:
            for impl in self.impl:
                atom.SetAtomMapNum(impl.NUCLEOPHILE_EAS_MAP_NUM)
                for template in templates:
                    try:
                        smarts = impl.generate(deepcopy(nucleophile_mol), template, atom, nucleophile)
                    except InvalidMolecule:
                        pass
                    else:
                        for rxn_str in smarts:
                            reactions.append(models.Reaction.from_mols(
                                impl.TYPE, rxn_str, template, nucleophile, atom.GetIdx()))

            atom.SetAtomMapNum(0)

        return reactions

    def _get_non_symmetric_atoms(self, mol):

        unique_pairs = map(set, zip(*mol.GetSubstructMatches(mol, uniquify=False)))
        Chem.Kekulize(mol)  # have to kekulize after substruct match...
        sorted_unique_pairs = map(sorted, map(list, unique_pairs))
        reduced_pairs = set(tuple(pair) for pair in sorted_unique_pairs)
        non_symmetric_atom_idxs = set(pair[0] for pair in reduced_pairs)

        # remove backbone atoms from non-symmetric atoms (this does nothing if mol is a sidechain)
        for backbone in self.backbones:
            temp_utils.clear_atom_map_nums(backbone)
            atom_idxs = set(chain.from_iterable(mol.GetSubstructMatches(backbone)))
            non_symmetric_atom_idxs = non_symmetric_atom_idxs - atom_idxs

        return non_symmetric_atom_idxs


class IntraMolecularReactionGenerator:
    STRING = 'intra_' + models.Reaction.STRING

    def __init__(self, impl=None):
        self.impl = impl or rxns.create_intramolecular_reactions()
        self.impl = temp_utils.to_list(self.impl)

    def generate(self, reacting_mol):
        reactions = []

        for impl in self.impl:
            try:
                for smarts in impl.generate(reacting_mol):
                    reactions.append(models.Reaction.from_mols(impl.TYPE, smarts, reacting_mol, None, None))
            except InvalidMolecule:
                pass

        return reactions


get_all_generator_strings = temp_utils.get_module_strings(__name__)

create_generator_from_string = temp_utils.create_factory_function_closure(__name__, 'generator')
