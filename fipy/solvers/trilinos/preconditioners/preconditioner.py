__docformat__ = 'restructuredtext'

__all__ = ["Preconditioner"]

class Preconditioner:
    """
    The base Preconditioner class.

    .. attention:: This class is abstract. Always create one of its subclasses.
    """

    def __init__(self):
        """
        Create a `Preconditioner` object.
        """
        if self.__class__ is Preconditioner:
            raise NotImplementedError("can't instantiate abstract base class")

    def _applyToSolver(self, solver, matrix):
        raise NotImplementedError
