class MissingMapNumber(Exception):
    """
    Raised when Mol is missing atom map numbers in SMILES string.
    """


class MergeError(Exception):
    """
    Raised when merging of two molecules fails.
    """
