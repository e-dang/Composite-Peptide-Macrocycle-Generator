from copy import deepcopy
from itertools import chain

from rdkit import Chem
from rdkit.Chem import AllChem

import cpmg.models as models
import cpmg.utils as temp_utils
import cpmg.repository as repo
import macrocycles.utils as utils
from cpmg.exceptions import InvalidMolecule
import cpmg.filters as filters


class SidechainModifier:

    def generate(self, args):
        """
        Method for creating new sidechain molecules from other sidechains by changing the methyl connection to a
        different connection type.
        """

        new_sidechains = []
        parent_sidechain, connections = args

        # replace the designated attachment position with each type of connection
        for connection in connections:
            for sidechain in Chem.ReplaceSubstructs(parent_sidechain.mol, models.SC_ATTACHMENT_POINT, connection.mol):
                new_sidechains.append(models.Sidechain.from_mol(sidechain, connection, parent_sidechain.shared_id))

        return new_sidechains


class MonomerGenerator:
    MAP_NUMS = (models.Backbone.MAP_NUM, models.Sidechain.MAP_NUM)

    def generate(self, args):
        """
        Method that takes a sidechain and backbones creates new monomers.
        """

        sidechain, backbones = args
        sidechain_mol = sidechain.mapped_mol

        temp_utils.clear_isotopes(sidechain_mol)

        monomers = []
        for backbone in backbones:
            monomer = utils.connect_mols(sidechain_mol, backbone.mol, map_nums=self.MAP_NUMS)
            monomers.append(models.Monomer.from_mol(monomer, backbone, sidechain))

        return monomers


class PeptideGenerator():

    MONOMER_NITROGEN_MAP_NUM = 1
    PEPTIDE_CARBON_MAP_NUM = 2
    MAP_NUMS = (MONOMER_NITROGEN_MAP_NUM, PEPTIDE_CARBON_MAP_NUM)

    # this reaction only works for alpha amino acids by design
    DECARBOXYLATE = AllChem.ReactionFromSmarts('[*:1]NC([*:2])C(=O)[OH]>>[*:1]NC([*:2])')

    def __init__(self, peptide_length):

        self.peptide_length = peptide_length

    def generate(self, monomers):
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
        if len(monomers) != self.peptide_length:  # c-term cap is present
            self.peptide = self.decarboxylate_c_term()
            has_c_cap = True

        return [models.Peptide.from_mol(self.peptide, has_c_cap, monomers)]

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
    # any primary amine or proline n-terminus, but no guanidine
    ELIGIBLE_NITROGENS = Chem.MolFromSmarts('[$([NH2]),$([NH;R]);!$([NH2]C(=O)*);!$([NH2]C(=[NH])[NH]*)]')
    PEPTIDE_NITROGEN_MAP_NUM = 2
    MAP_NUMS = (models.Template.OLIGOMERIZATION_MAP_NUM, PEPTIDE_NITROGEN_MAP_NUM)

    def generate(self, args):
        """
        Method that takes a peptide molecule and combines it with each type of template molecule to form
        template-peptide oligomers.
        """

        peptide, templates = args

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


class InterMolecularReactionGenerator:

    def __init__(self, impl):
        self.impl = impl
        self.backbones = list(map(lambda x: x.mol, repo.create_backbone_repository().load()))

    @filters.pka_filter
    @filters.regiosqm_filter
    def generate(self, args):
        nucleophile, templates = args
        nucleophile_mol = nucleophile.mol

        non_symmetric_atom_idxs = self._get_non_symmetric_atoms(nucleophile_mol)

        reactions = []
        for atom in [nucleophile_mol.GetAtomWithIdx(atom_idx) for atom_idx in non_symmetric_atom_idxs]:
            atom.SetAtomMapNum(self.impl.NUCLEOPHILE_EAS_MAP_NUM)
            for template in templates:
                try:
                    smarts = self.impl.generate(deepcopy(nucleophile_mol), template, atom, nucleophile)
                except InvalidMolecule:
                    pass
                else:
                    for rxn_str in smarts:
                        reactions.append(models.Reaction.from_mols(
                            self.impl.TYPE, rxn_str, template, nucleophile, atom.GetIdx()))

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
