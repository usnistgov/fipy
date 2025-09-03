__docformat__ = 'restructuredtext'

__all__ = ["JacobiPreconditioner"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

from scipy.sparse import diags
from scipy.sparse.linalg import LinearOperator, spsolve

from .scipyPreconditioner import ScipyPreconditioner

class JacobiPreconditioner(ScipyPreconditioner):
    """Jacobi preconditioner for :class:`~fipy.solvers.scipy.scipySolver.ScipySolver`.

    Wrapper class for :func:`scipy.sparse.linalg.spsolve` with `matrix`
    diagonal.
    Adapted from https://stackoverflow.com/q/46876951/2019542.
    """
    
    def _applyToMatrix(self, matrix):
        P = diags(matrix.diagonal()).tocsc()
        Mx = lambda x: spsolve(P, x)

        return LinearOperator(matrix.shape, Mx), matrix
