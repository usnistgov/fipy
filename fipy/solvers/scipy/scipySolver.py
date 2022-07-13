from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

__all__ = []

from scipy.sparse import linalg

from fipy.matrices.scipyMatrix import _ScipyMeshMatrix
from fipy.solvers.solver import Solver
from fipy.tools import numerix

class _ScipySolver(Solver):
    """
    The base `ScipySolver` class.

    .. attention:: This class is abstract. Always create one of its subclasses.
    """

    def __init__(self, tolerance=1e-10, criterion="default",
                 iterations=1000, precon=None):
        """
        Create a `Solver` object.

        Parameters
        ----------
        tolerance : float
            Required error tolerance.
        criterion : {'default', 'unscaled', 'RHS', 'matrix', 'initial'}
            Interpretation of ``tolerance``.
            See :ref:`CONVERGENCE` for more information.
        iterations : int
            Maximum number of iterative steps to perform.
        precon
            Preconditioner to use.  Not all solver suites support
            preconditioners.
        """
        super(_ScipySolver, self).__init__(tolerance=tolerance, criterion=criterion,
                                           iterations=iterations, precon=precon)

    @property
    def _matrixClass(self):
        return _ScipyMeshMatrix

    def _rhsNorm(self, L, x, b):
        return numerix.L2norm(b)

    def _matrixNorm(self, L, x, b):
        return linalg.norm(L.matrix, ord=numerix.inf)

    def _residualVectorAndNorm(self, L, x, b):
        residualVector = L * x - b

        return residualVector, numerix.L2norm(residualVector)

    def _solve(self):

         if self.var.mesh.communicator.Nproc > 1:
             raise Exception("SciPy solvers cannot be used with multiple processors")

         self.var[:] = numerix.reshape(self._solve_(self.matrix, self.var.ravel(), numerix.array(self.RHSvector)), self.var.shape)
