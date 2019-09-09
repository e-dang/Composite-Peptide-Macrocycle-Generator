"""
Written by Eric Dang.
github: https://github.com/e-dang
email: edang830@gmail.com
"""


class MissingMapNumberError(Exception):
    """
    Raised when Mol is missing atom map numbers in SMILES string.
    """


class AtomSearchError(Exception):
    """
    Raised when unable to find atom that meets specified criteria.
    """


class InvalidSmilesString(Exception):
    """
    Raised when unable to convert a SMILES/SMARTS string to rdkit Mol.
    """


class DataNotFoundError(Exception):
    """
    Raised when the MongoDataBase returns an empty result.
    """


class MergeError(Exception):
    """
    Raised when merging of two molecules fails.
    """


class FailedEmbeddingError(Exception):
    """
    Raise when embedding a molecule for conformation generation using ETKDGv2 parameters and random coordinates both
    fail to produce a starting conformation.
    """


class CTermNotFound(Exception):
    """
    Raised when searching for the c-terminus carboxyl group on a macrocycle fails.
    """
