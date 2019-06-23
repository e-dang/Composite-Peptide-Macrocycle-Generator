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


class WritingJsonError(Exception):
    """
    Raised when an error has occured while writing data to a json file.
    """


class WritingTxtError(Exception):
    """
    Raised when an error has occured while writing data to a txt file.
    """


class SavingMongoError(Exception):
    """
    Raised when an error has occured while saving data to the Mongo database.
    """


class SavingSQLError(Exception):
    """
    Raised when an error has occured while saving data to the SQL database.
    """
