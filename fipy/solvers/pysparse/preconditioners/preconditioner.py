__all__ = ["Preconditioner"]

class Preconditioner:
    """
    Base preconditioner class

    .. attention:: This class is abstract. Always create one of its subclasses.
    """

    def __init__(self):
        """
        Create a `Preconditioner` object.
        """
        if self.__class__ is Preconditioner:
            raise NotImplementedError("can't instantiate abstract base class")

    def _applyToMatrix(self, matrix):
        """
        Returns the function used for PySparse
        preconditioning.
        """
        raise NotImplementedError
