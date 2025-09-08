__docformat__ = 'restructuredtext'

from fipy.solvers.preconditioner import MatrixModifyingPreconditioner

__all__ = ["ScipyPreconditioner"]

class ScipyPreconditioner(MatrixModifyingPreconditioner):
    """Base class for preconditioners for :class:`~fipy.solvers.scipy.scipySolver.ScipySolver`.

    .. attention:: This class is abstract. Always create one of its subclasses.
    """

    pass
