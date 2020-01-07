from rdkit import Chem
from abc import ABC, abstractmethod


class IMolecule(ABC):
    """
    Interface for classes that represent molecules. These types of molecules are ones that need to be tagged with atom
    map numbers or have certain modifications at certain positions that are hard to pinpoint with code and are more
    easily done by manually doing it in a third party software application such as ChemDraw.
    """

    @property
    @abstractmethod
    def name(self):
        """
        Abstract property that returns the name of the molecule.
        """

    @property
    @abstractmethod
    def type(self):
        """
        Abstract property that returns the type of the molecule (used for organizing data in this project).
        """

    @property
    @abstractmethod
    def kekule(self):
        """
        Property that returns the kekulized SMILES string of the molecule associated with the specific derived class.
        """

    @property
    def mol(self):
        """
        Property that returns the RDKit Mol associated with the specific derived class.
        """
        return Chem.MolFromSmiles(self.kekule)


class ITemplateMol(IMolecule):
    """
    Interface for template molecules.
    """

    OLIGOMERIZATION_MAP_NUM = 1  # for tagging the carbon on the aldehyde where the peptide will oligomerize to template
    EAS_CARBON_MAP_NUM = 2  # for tagging the carbon atom that will close the macrocycle ring in EAS reaction
    PS_CARBON_ALDEHYDE_MAP_NUM = 7  # for tagging the carbon on the aldehye required for pictet spangler reaction
    PS_OXYGEN_ALDEHYDE_MAP_NUM = 8  # for tagging the oxygen on the aldehyde requried for pictet spangler reaction
    TEMPLATE_PS_REACTING_ATOM_MAP_NUM = 9  # for tagging the reacting atom in template pictet spangler reaction
    TEMPLATE_PS_NITROGEN_MAP_NUM = 10  # for tagging the peptide nitrogen in template pictet spangler reaction
    REACTION_WILDCARD1_MAP_NUM = 50  # for a wildcard atom
    REACTION_WILDCARD2_MAP_NUM = 51  # for a wildcard atom

    @property
    def type(self):
        return 'template'

    @property
    @abstractmethod
    def oligomerization_kekule(self):
        """
        Abstract property that returns the atom mapped kekulized SMILES string of the template molecule that has been
        edited to contain an atom map number equal to OLIGOMERIZATION_MAP_NUM at the aldehyde carbon that is to be
        connected to the N-terminus nitrogen of a peptide to form a template-peptide oligomer. Other changes to the
        original template kekule SMILES string may have been made for convenience as well.
        """

    @property
    def oligomerization_mol(self):
        """
        Generates an RDKit Mol from the oligomerization_kekule SMILES string.

        Returns:
            RDKit Mol: The RDKit Mol object.
        """
        return Chem.MolFromSmiles(self.oligomerization_kekule)


class CinnamoylTemplate1(ITemplateMol):
    """
    Generation 1 of the cinnamoyl template molecule.
    """

    @property
    def name(self):
        return 'temp1'

    @property
    def kekule(self):
        return 'CC(C)(C)OC(=O)OC/C=C/C1=CC(CCC(=O)ON2C(=O)CCC2=O)=CC=C1'

    @property
    def friedel_crafts_kekule(self):
        return f'[*:{self.REACTION_WILDCARD1_MAP_NUM}]/C=C/[CH3:{self.EAS_CARBON_MAP_NUM}]'

    @property
    def tsuji_trost_kekule(self):
        return self.friedel_crafts_kekule

    @property
    def pyrrolo_indolene_kekule(self):
        return self.friedel_crafts_kekule

    @property
    def oligomerization_kekule(self):
        return f'C/C=C/C1=CC(CC[CH:{self.OLIGOMERIZATION_MAP_NUM}]=O)=CC=C1'


class CinnamoylTemplate2(ITemplateMol):
    """
    Generation 2 of the cinnamoyl template molecule.
    """

    @property
    def name(self):
        return 'temp2'

    @property
    def kekule(self):
        return 'CC(C)(C)OC(=O)OC/C=C/C1=CC(C[C@@H](CC=O)C(=O)ON2C(=O)CCC2=O)=C(F)C=C1'

    @property
    def pictet_spangler_kekule(self):
        return f'[O:{self.PS_OXYGEN_ALDEHYDE_MAP_NUM}]=[CH1:{self.PS_CARBON_ALDEHYDE_MAP_NUM}]C[C@H](C[*:{self.REACTION_WILDCARD1_MAP_NUM}])[CH1:{self.OLIGOMERIZATION_MAP_NUM}]=O'

    @property
    def oligomerization_kekule(self):
        return f'C/C=C/C1=CC=C(F)C(C[C@@H](CC=O)[CH:{self.OLIGOMERIZATION_MAP_NUM}]=O)=C1'


class CinnamoylTemplate3(ITemplateMol):
    """
    Generation 3 of the cinnamoyl template molecule.
    """

    @property
    def name(self):
        return 'temp3'

    @property
    def kekule(self):
        return 'C#CCCC[C@@](Cc1cc(/C=C/COC(OC(C)(C)C)=O)ccc1)(C=O)CC(ON2C(CCC2=O)=O)=O'

    @property
    def pictet_spangler_kekule(self):
        return f'C#CCCC[C@](C[*:{self.REACTION_WILDCARD1_MAP_NUM}])(C[CH1:{self.OLIGOMERIZATION_MAP_NUM}]=O)[CH1:{self.PS_CARBON_ALDEHYDE_MAP_NUM}]=[O:{self.PS_OXYGEN_ALDEHYDE_MAP_NUM}]'

    @property
    def template_pictet_spangler_kekule(self):
        return f'C#CCCC[C@@]([CH1:{self.PS_CARBON_ALDEHYDE_MAP_NUM}]=[O:{self.PS_OXYGEN_ALDEHYDE_MAP_NUM}])(CC(=O)[NH1:{self.TEMPLATE_PS_NITROGEN_MAP_NUM}][*:{self.REACTION_WILDCARD1_MAP_NUM}])CC1=[CH1:{self.TEMPLATE_PS_REACTING_ATOM_MAP_NUM}]C=CC([*:{self.REACTION_WILDCARD2_MAP_NUM}])=C1'

    @property
    def oligomerization_kekule(self):
        return f'C#CCCC[C@](C=O)(CC1=CC=CC(/C=C/C)=C1)C[CH:{self.OLIGOMERIZATION_MAP_NUM}]=O'


class IConnectionMol(IMolecule):
    """
    Interface for molecules to be used as a connection between a heterocycle and a backbone in order to form a monomer.
    """

    @property
    def type(self):
        return 'connection'


class EthylConnection(IConnectionMol):
    """
    A connection molecule that is two carbon atoms long.
    """

    @property
    def name(self):
        return 'ethyl'

    @property
    def kekule(self):
        return 'CC'


class IBackBoneMol(IMolecule):
    """
    Interface for molecules that are used as backbone structures for forming monomers. These molecules should have the
    position on their backbone to which a sidechain will be attached marked with the atom map number
    OLIGOMERIZATION_MAP_NUM.
    """

    OLIGOMERIZATION_MAP_NUM = 1

    @property
    def type(self):
        return 'backbone'

    @property
    @abstractmethod
    def tagged_kekule(self):
        """
        Abstract property that returns the atom mapped kekulized SMILES string of the backbone molecule needed for
        monomer generation.
        """

    @property
    def tagged_mol(self):
        """
        Generates an RDKit Mol from the tagged_kekule SMILES string.

        Returns:
            RDKit Mol: The RDKit Mol object.
        """
        return Chem.MolFromSmiles(self.tagged_kekule)


class AlphaBackBone(IBackBoneMol):
    """
    The backbone structure of a alpha amino acid.
    """

    @property
    def name(self):
        return 'alpha'

    @property
    def kekule(self):
        return 'NCC(=O)[OH]'

    @property
    def tagged_kekule(self):
        return f'N[CH2:{self.OLIGOMERIZATION_MAP_NUM}]C(=O)O'


class Beta2BackBone(IBackBoneMol):
    """
    The backbone structure of a beta_2 amino acid, where the connection point is adjacent to the nitrogen.
    """

    @property
    def name(self):
        return 'beta2'

    @property
    def kekule(self):
        return 'NCCC(=O)[OH]'

    @property
    def tagged_kekule(self):
        return f'N[CH2:{self.OLIGOMERIZATION_MAP_NUM}]CC(=O)O'


class Beta3BackBone(IBackBoneMol):
    """
    The backbone structure of a beta_3 amino acid, where the connection point is adjacent to the carboxyl group.
    """

    @property
    def name(self):
        return 'beta3'

    @property
    def kekule(self):
        return 'NCCC(=O)[OH]'

    @property
    def tagged_kekule(self):
        return f'NC[CH2:{self.OLIGOMERIZATION_MAP_NUM}]C(=O)O'
