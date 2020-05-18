class MergeError(Exception):
    """
    Raised when merging of two molecules fails.
    """


class InvalidMolecule(Exception):
    """
    Raised when an operation on a molecule cannot be performed due to the molecule missing a feature required for that
    operation.
    """


class InvalidPrediction(Exception):
    """
    Raised when importing a prediction that fails the sanity test.
    """


class ParallelismAlreadySet(Exception):
    """
    Raised when trying to change the parallelism level when it has already been set.
    """


class TimerAlreadyStarted(Exception):
    """
    Raised when start() has been called on timer class when it has already been started.
    """


class TimerNotStarted(Exception):
    """
    Raised when calling a function on Timer class that determines how much time has been elapsed when the Timer hasn't been started yet.
    """


class InvalidChunkSize(Exception):
    """
    Raised when a negative or 0 chunk size is specified
    """
