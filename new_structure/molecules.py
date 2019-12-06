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

    _OLIGOMERIZATION_MAP_NUM = 1

    # @property
    # @abstractmethod
    # def reaction_kekule(self):
    #     """
    #     Abstract property that returns an edited atom mapped version of the kekulized SMILES string of the template
    #     molecule. The required edits made to the kekule SMILES string are that
    #     """

    @property
    @abstractmethod
    def oligomerization_kekule(self):
        """
        Abstract property that returns the atom mapped kekulized SMILES string of the template molecule that has been
        edited to contain an atom map number equal to _OLIGOMERIZATION_MAP_NUM at the aldehyde carbon that is to be
        connected to the N-terminus nitrogen of a peptide to form a template-peptide oligomer. Other changes to the
        original template kekule SMILES string may have been made for convenience as well.
        """


class CinnamoylTemplate1(ITemplateMol):
    """
    Generation 1 of the cinnamoyl template molecule.
    """

    @property
    def name(self):
        return 'temp1'

    @property
    def type(self):
        return 'template'

    @property
    def kekule(self):
        return 'CC(C)(C)OC(=O)OC/C=C/C1=CC(CCC(=O)ON2C(=O)CCC2=O)=CC=C1'

    # @property
    # def reaction_kekule(self):
    #     return 'C(=C/[CH3:2])\C1=CC=CC(CC[*:1])=C1'

    @property
    def oligomerization_kekule(self):
        return f'C/C=C/C1=CC(CC[CH:{self._OLIGOMERIZATION_MAP_NUM}]=O)=CC=C1'


class CinnamoylTemplate2(ITemplateMol):
    """
    Generation 2 of the cinnamoyl template molecule.
    """

    @property
    def name(self):
        return 'temp2'

    @property
    def type(self):
        return 'template'

    @property
    def kekule(self):
        return 'CC(C)(C)OC(=O)OC/C=C/C1=CC(CCC(=O)ON2C(=O)CCC2=O)=CC=C1'

    # @property
    # def reaction_kekule(self):
    #     return 'O=CC[C@H](CC1=CC(/C=C/[CH3:2])=CC=C1F)[*:1]'

    @property
    def oligomerization_kekule(self):
        return f'C/C=C/C1=CC=C(F)C(C[C@@H](CC=O)[CH:{self._OLIGOMERIZATION_MAP_NUM}]=O)=C1'


class CinnamoylTemplate3(ITemplateMol):
    """
    Generation 3 of the cinnamoyl template molecule.
    """

    @property
    def name(self):
        return 'temp3'

    @property
    def type(self):
        return 'template'

    @property
    def kekule(self):
        return 'CC(C)(C)OC(=O)OC/C=C/C1=CC(CCC(=O)ON2C(=O)CCC2=O)=CC=C1'

    # @property
    # def reaction_kekule(self):
    #     return 'C#CCCC[C@](C=O)(CC1=CC=CC(/C=C/[CH3:2])=C1)C[*:1]'

    @property
    def oligomerization_kekule(self):
        return f'C#CCCC[C@](C=O)(CC1=CC=CC(/C=C/C)=C1)C[CH:{self._OLIGOMERIZATION_MAP_NUM}]=O'


TEMPLATE_MOLS = [CinnamoylTemplate1(), CinnamoylTemplate2(), CinnamoylTemplate3()]


class IConnectionMol(IMolecule):
    """
    Interface for molecules to be used as a connection between a heterocycle and a backbone in order to form a monomer.
    """

    _OLIGOMERIZATION_MAP_NUM = 1

    @property
    def type(self):
        return 'connection'

    @property
    @abstractmethod
    def tagged_kekule(self):
        """
        Abstract property that returns the atom mapped kekulized SMILES string of the connection molecule needed for
        sidechain generation.
        """


class MethylConnection(IConnectionMol):
    """
    A connection molecule that is a single carbon atom long.
    """

    @property
    def name(self):
        return 'methyl'

    @property
    def kekule(self):
        return 'C'

    @property
    def tagged_kekule(self):
        return f'[CH4:{self._OLIGOMERIZATION_MAP_NUM}]'

    @property
    def mol(self):
        return Chem.MolFromSmiles(self.tagged_kekule)


class EthylConnection(IConnectionMol):
    """
    A connection molecule that is two carbon atoms long.
    """

    @property
    def name(self):
        return 'methyl'

    @property
    def kekule(self):
        return 'CC'

    @property
    def tagged_kekule(self):
        return f'C[CH3:{self._OLIGOMERIZATION_MAP_NUM}]'

    @property
    def mol(self):
        return Chem.MolFromSmiles(self.tagged_kekule)


CONNECTION_MOLS = [MethylConnection(), EthylConnection()]


class IBackBoneMol(IMolecule):
    """
    Interface for molecules that are used as backbone structures for forming monomers. These molecules should have the
    position on their backbone to which a sidechain will be attached marked with the atom map number
    _OLIGOMERIZATION_MAP_NUM.
    """

    _OLIGOMERIZATION_MAP_NUM = 1

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


class AlphaBackBone(IBackBoneMol):
    """
    The backbone structure of a alpha amino acid.
    """

    @property
    def name(self):
        return 'alpha'

    @property
    def kekule(self):
        return 'NCC(=O)O'

    @property
    def tagged_kekule(self):
        return f'N[CH2:{self._OLIGOMERIZATION_MAP_NUM}]C(=O)O'


class Beta2BackBone(IBackBoneMol):
    """
    The backbone structure of a beta_2 amino acid, where the connection point is adjacent to the nitrogen.
    """

    @property
    def name(self):
        return 'alpha'

    @property
    def kekule(self):
        return 'NCCC(=O)O'

    @property
    def tagged_kekule(self):
        return f'N[CH2:{self._OLIGOMERIZATION_MAP_NUM}]CC(=O)O'


class Beta3BackBone(IBackBoneMol):
    """
    The backbone structure of a beta_3 amino acid, where the connection point is adjacent to the carboxyl group.
    """

    @property
    def name(self):
        return 'alpha'

    @property
    def kekule(self):
        return 'NCCC(=O)O'

    @property
    def tagged_kekule(self):
        return f'NC[CH2:{self._OLIGOMERIZATION_MAP_NUM}]C(=O)O'


BACKBONE_MOLS = [AlphaBackBone(), Beta2BackBone(), Beta3BackBone()]
