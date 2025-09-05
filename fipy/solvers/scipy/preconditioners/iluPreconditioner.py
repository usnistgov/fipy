__docformat__ = 'restructuredtext'

__all__ = ["ILUPreconditioner"]

from scipy.sparse.linalg import LinearOperator, spilu

from .scipyPreconditioner import ScipyPreconditioner

class ILUPreconditioner(ScipyPreconditioner):
    """Incomplete LU preconditioner for :class:`~fipy.solvers.scipy.scipySolver.ScipySolver`.

    Wrapper class for :func:`scipy.sparse.linalg.spilu`.
    Adapted from https://stackoverflow.com/q/46876951/2019542.
    """
    
    def _applyToMatrix(self, matrix):
        ilu = spilu(matrix.tocsc())
        Mx = lambda x: ilu.solve(x)

        return LinearOperator(matrix.shape, Mx), matrix
