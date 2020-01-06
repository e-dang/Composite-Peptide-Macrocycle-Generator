import exceptions
from abc import ABC, abstractmethod
from copy import deepcopy
from itertools import chain

from rdkit import Chem
from rdkit.Chem import AllChem

import decorators
import molecules
import reactions
import utils


class IGenerator(ABC):
    """
    Interface for classes that operate on molecular structure(s) in order to generate new structures.
    """

    @abstractmethod
    def generate(self, args):
        """
        Abstract method for applying the specific molecular transformation handled by the derived class on the molecules
        passed in as arguments.

        Args:
            args (iterable): An iterable containing the molecules and associated data that are to be operated on.
        """


class SideChainConnectionModifier(IGenerator):
    """
    Implementation of an IGenerator that takes a sidechain with a methyl connection and creates a set of sidechains by
    replacing the methyl connection with all other connection types defined in molecules.py. The new sidechains are
    assigned a new id by the id_iterator used to initialize this generator.
    """

    def __init__(self, id_iterator):
        self.id_iterator = id_iterator

    def generate(self, args):

        self.new_sidechains = {}
        self.parent_sidechain = args
        parent_sidechain = Chem.Mol(self.parent_sidechain['binary'])

        # methyls with carbon 13 are marked as pure methyl groups not sites of attachment
        patt = Chem.MolFromSmarts('[CH3;!13CH3]')

        # replace the designated attachment position with each type of connection
        for connection in utils.get_connections():
            for sidechain in Chem.ReplaceSubstructs(parent_sidechain, patt, connection.mol):
                Chem.SanitizeMol(sidechain)
                binary = sidechain.ToBinary()
                Chem.Kekulize(sidechain)
                self.new_sidechains[Chem.MolToSmiles(sidechain, kekuleSmiles=True)] = (binary, connection.name)

        return self.format_data()

    def format_data(self):
        """
        Helper method that fills in a dict object for each new sidechain with the necessary data associated with that
        sidechain needed for record keeping.

        Returns:
            list: A list containing the newly created sidechain dicts.
        """

        data = []
        for kekule, (binary, connection) in self.new_sidechains.items():
            data.append({'_id': self.id_iterator.get_next(),
                         'type': 'sidechain',
                         'binary': binary,
                         'kekule': kekule,
                         'connection': connection,
                         'shared_id': self.parent_sidechain['_id']})

        return data


class MonomerGenerator(IGenerator):
    """
    Implementation of an IGenerator that takes a sidechain and backbone molecule and creates a monomer
    by attaching the alkyl portion of the sidechain to the designated position on the backbone molecule.
    """

    BB_MAP_NUM = molecules.IBackBoneMol.OLIGOMERIZATION_MAP_NUM
    SC_MAP_NUM = 2
    MAP_NUMS = (BB_MAP_NUM, SC_MAP_NUM)

    def __init__(self, index_iterator):
        self.index_iterator = index_iterator

    def generate(self, args):

        self.monomers = {}
        self.sidechain = args
        sidechain = Chem.Mol(self.sidechain['binary'])
        required = bool(AllChem.CalcNumAromaticRings(sidechain))

        # find candidate attachment point at terminal end of alkyl chains
        matches = sidechain.GetSubstructMatches(Chem.MolFromSmarts('[CH3;!13CH3]'))
        if len(matches) != 1:
            raise exceptions.InvalidMolecule(
                f'Sidechain {Chem.MolToSmiles(sidechain)} is missing an attachment point')

        # set atom map number for attachment point of sidechain
        for atom_idx in chain.from_iterable(matches):
            sidechain.GetAtomWithIdx(atom_idx).SetAtomMapNum(self.SC_MAP_NUM)

        # change isotopes of non-candidate attachment points back to normal isotope
        for atom_idx in chain.from_iterable(sidechain.GetSubstructMatches(Chem.MolFromSmarts('[13CH3]'))):
            sidechain.GetAtomWithIdx(atom_idx).SetIsotope(12)

        # connect monomer and backbone
        for backbone in utils.get_backbones():
            monomer = utils.connect_mols(deepcopy(sidechain), backbone.tagged_mol, map_nums=self.MAP_NUMS)
            binary = monomer.ToBinary()
            Chem.Kekulize(monomer)
            self.monomers[Chem.MolToSmiles(monomer, kekuleSmiles=True)] = (binary, required, backbone.name)

        return self.format_data()

    def format_data(self):
        """
        Helper method that fills in a dict object for each new monomer with the necessary data associated with that
        monomer needed for record keeping.

        Returns:
            list: A list containing the newly created monomer dicts.
        """

        data = []
        for i, (kekule, (binary, required, backbone)) in enumerate(self.monomers.items()):
            data.append({'_id': self.sidechain['_id'].upper() + str(i),
                         'type': 'monomer',
                         'binary': binary,
                         'kekule': kekule,
                         'index': self.index_iterator.get_next(),
                         'required': required,
                         'backbone': backbone,
                         'sidechain': self.sidechain['shared_id']})

        return data


class PeptideGenerator(IGenerator):

    BACKBONES = utils.get_hashed_backbones()
    MONOMER_NITROGEN_MAP_NUM = 1
    PEPTIDE_CARBON_MAP_NUM = 2
    MAP_NUMS = (MONOMER_NITROGEN_MAP_NUM, PEPTIDE_CARBON_MAP_NUM)

    def generate(self, args):

        self.monomers = args

        # begin conneting each monomer in monomers
        for i, monomer_doc in enumerate(self.monomers):

            self.monomer = Chem.Mol(monomer_doc['binary'])
            self.backbone = self.BACKBONES[monomer_doc['backbone']]

            # start peptide with first monomer
            if i == 0:
                self.pep_size = 1
                self.peptide = self.monomer
                self.backbone_prev = self.backbone
                continue

            # assign atom map numbers
            monomer_old_attach = self.tag_monomer_n_term()
            carboxyl_atom, pep_old_attach = self.tag_peptide_c_term()

            # remove oxygen atom from carboxyl
            self.peptide = Chem.RWMol(self.peptide)
            self.peptide.RemoveAtom(carboxyl_atom)

            # connect peptide and monomer
            self.peptide = utils.connect_mols(self.peptide, self.monomer, map_nums=self.MAP_NUMS)
            self.pep_size += 1
            pep_old_attach.SetAtomMapNum(0)
            monomer_old_attach.SetAtomMapNum(0)
            self.backbone_prev = self.backbone

        return self.format_data()

    def tag_monomer_n_term(self):

        for atom_idx in self.monomer.GetSubstructMatch(self.backbone):
            atom = self.monomer.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() == 'N' and atom.GetTotalNumHs() != 0:
                atom.SetAtomMapNum(self.MONOMER_NITROGEN_MAP_NUM)
                return atom

    def tag_peptide_c_term(self):

        flag = False
        for pair in self.peptide.GetSubstructMatches(self.backbone_prev):
            for atom_idx in pair:
                atom = self.peptide.GetAtomWithIdx(atom_idx)
                if atom.GetSymbol() == 'O' and atom.GetTotalNumHs() == 1:  # finds protonated oxygen on carboxyl group
                    carboxyl_atom = atom_idx
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

        return carboxyl_atom, attachment_atom

    def format_data(self):
        monomer_data = [{key: value for key, value in monomer.items() if key in ('_id', 'sidechain')}
                        for monomer in self.monomers]
        pep_id = ''.join([monomer['_id'] for monomer in monomer_data])

        binary = self.peptide.ToBinary()
        Chem.Kekulize(self.peptide)
        return [{'_id': pep_id,
                 'type': 'peptide',
                 'binary': binary,
                 'kekule': Chem.MolToSmiles(self.peptide, kekuleSmiles=True),
                 'monomers': monomer_data}]


class TemplatePeptideGenerator(IGenerator):

    # any primary amine or proline n-terminus
    ELIGIBLE_NITROGENS = Chem.MolFromSmarts('[$([NH2]),$([NH;R]);!$([NH2]C(=O)*)]')
    TEMPLATE_CARBON_MAP_NUM = molecules.ITemplateMol.OLIGOMERIZATION_MAP_NUM
    PEPTIDE_NITROGEN_MAP_NUM = 2
    MAP_NUMS = (TEMPLATE_CARBON_MAP_NUM, PEPTIDE_NITROGEN_MAP_NUM)

    def generate(self, args):

        self.template_peptides = {}
        self.peptide = args

        # for each template and eligible nitrogen in peptide form a connection
        for template in utils.get_templates():
            peptide = Chem.Mol(self.peptide['binary'])
            for atom_idx in chain.from_iterable(peptide.GetSubstructMatches(self.ELIGIBLE_NITROGENS)):
                atom = peptide.GetAtomWithIdx(atom_idx)
                atom.SetAtomMapNum(self.PEPTIDE_NITROGEN_MAP_NUM)
                template_peptide = utils.connect_mols(template.oligomerization_mol, peptide, map_nums=self.MAP_NUMS)
                atom.SetAtomMapNum(0)

                binary = template_peptide.ToBinary()
                Chem.Kekulize(template_peptide)
                self.template_peptides[Chem.MolToSmiles(template_peptide, kekuleSmiles=True)] = (binary,
                                                                                                 template.name)
        return self.format_data()

    def format_data(self):

        data = []
        for i, (kekule, (binary, template)) in enumerate(self.template_peptides.items()):
            data.append({'_id': self.peptide['_id'] + str(i),
                         'type': 'template_peptide',
                         'binary': binary,
                         'kekule': kekule,
                         'template': template,
                         'peptide': {'_id': self.peptide['_id'],
                                     'monomers': self.peptide['monomers']}})

        return data


class MacrocycleGenerator(IGenerator):


    @decorators.apply_stereochemistry
    @decorators.methylate
    @decorators.carboxyl_to_amide
    def generate(self, args):

        self.reactant, self.reactions = args
        reactants = [Chem.Mol(self.reactant['binary'])]
        rxns = [AllChem.ChemicalReaction(reaction['binary']) for reaction in self.reactions]
        # rxns = [(AllChem.ChemicalReaction(reaction['binary']), reaction['type']) for reaction in self.reactions]

        for rxn in rxns:
            self.macrocycles = {}
            for reactant in reactants:
                for macrocycle in chain.from_iterable(rxn.RunReactants((reactant,))):
                    try:
                        Chem.SanitizeMol(macrocycle)
                    except ValueError:
                        # most likely there are 2 sidechains where one contains the other as a substructure which
                        # causes this exception to be raised. This is expect to occur quite often.
                        continue
                    for atom in chain.from_iterable(rxn.GetReactingAtoms()):
                        atom = macrocycle.GetAtomWithIdx(atom)
                        if atom.GetIsAromatic():
                            atom.SetProp('_protected', '1')
                    self.macrocycles[macrocycle.ToBinary()] = macrocycle
            # if rxn_type in ('pictet_spangler', 'template_pictet_spangler', 'unmasked_aldehyde_cyclization') and len(list(self.macrocycles.values())) == 0:
            #     continue
            reactants = list(self.macrocycles.values())

        self.macrocycles = {}  # final results at this point are stored in reactants variable
        for macrocycle in reactants:
            binary = macrocycle.ToBinary()
            Chem.Kekulize(macrocycle)
            self.macrocycles[Chem.MolToSmiles(macrocycle, kekuleSmiles=True)] = binary

        return self.format_data()

    def format_data(self):
        rxns = []
        sc_id, rxn_idx = '', ''
        for reaction in self.reactions:
            rxns.append(reaction['_id'])
            try:
                rxn_idx += str(reaction['rxn_atom_idx'])
                sc_id += reaction['sidechain']
            except KeyError:
                pass

        data = []
        for i, (kekule, binary) in enumerate(self.macrocycles.items()):
            data.append({'_id': self.reactant['_id'] + sc_id + rxn_idx + str(i),
                         'type': 'macrocycle',
                         'binary': binary,
                         'kekule': kekule,
                         'has_confs': False,
                         'modifications': '',
                         'template': self.reactant['template'],
                         'template_peptide': self.reactant['_id'],
                         'reactions': rxns})
        return data


class MacrocycleConformerGenerator(IGenerator):

    def generate(self, args):
        pass


class UniMolecularReactionGenerator(IGenerator):

    def generate(self, args):

        self.reacting_mol, self.reaction = args
        self.reaction(self.reacting_mol)
        if self.reaction:
            self.reaction.create_reaction()
            return self.format_data()

        return []

    def format_data(self):

        data = [{'_id': self.reacting_mol.name + self.reaction.name[:2],
                 'type': self.reaction.name,
                 'binary': self.reaction.binary,
                 'smarts': self.reaction.smarts,
                 'template': self.reacting_mol.name}]

        return data


class BiMolecularReactionGenerator(IGenerator):

    REACTING_MOL_EAS_MAP_NUM = reactions.IReaction.REACTING_MOL_EAS_MAP_NUM

    @decorators.pka_filter
    @decorators.regiosqm_filter
    def generate(self, args):

        self.reactions = {}
        self.reacting_mol, self.reaction = args
        reacting_mol = Chem.MolFromSmiles(self.reacting_mol['kekule']) # need to do this since filter predictions are
                                                                       # based on atom indices and binary doesn't
                                                                       # maintain the same order

        # get non-symmetric atoms
        unique_pairs = map(set, zip(*reacting_mol.GetSubstructMatches(reacting_mol, uniquify=False)))
        Chem.Kekulize(reacting_mol) # have to kekulize after substruct match...
        sorted_unique_pairs = map(sorted, map(list, unique_pairs))
        reduced_pairs = set(tuple(pair) for pair in sorted_unique_pairs)
        non_symmetric_atom_idxs = [reacting_mol.GetAtomWithIdx(pair[0]) for pair in reduced_pairs]

        # create reactions
        for atom in non_symmetric_atom_idxs:
            atom.SetAtomMapNum(self.REACTING_MOL_EAS_MAP_NUM)
            for template in utils.get_templates():
                self.reaction(deepcopy(reacting_mol), template, atom)
                if self.reaction:
                    self.reaction.create_reaction()
                    self.reactions[self.reaction.smarts] = (self.reaction.binary, atom.GetIdx(),
                                                            self.reaction.name, self.reaction.applicable_template)

            atom.SetAtomMapNum(0)

        return self.format_data()

    def format_data(self):

        data = []
        for i, (smarts, (binary, rxn_atom_idx, rxn_type, template)) in enumerate(self.reactions.items()):
            reacting_mol_id = self.reacting_mol['shared_id'] if self.reacting_mol['type'] == 'sidechain' \
                                else self.reacting_mol['_id']
            data.append({'_id': self.reacting_mol['_id'] + str(i) + rxn_type[:2],
                         'type': rxn_type,
                         'binary': binary,
                         'smarts': smarts,
                         'rxn_atom_idx': rxn_atom_idx,
                         'template': template,
                         'reacting_mol': reacting_mol_id})

        return data
