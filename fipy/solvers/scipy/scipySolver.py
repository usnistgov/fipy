from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

__all__ = []

from fipy.matrices.scipyMatrix import _ScipyMeshMatrix
from fipy.solvers.solver import Solver
from fipy.tools import numerix

class _ScipySolver(Solver):
    """
    The base `ScipySolver` class.

    .. attention:: This class is abstract. Always create one of its subclasses.
    """

    criteria = {
        "default": None,
        "RHS": None
    }

    def __init__(self, tolerance=1e-10, criterion="initial",
                 iterations=1000, precon=None):
        """
        Create a `Solver` object.

        Parameters
        ----------
        tolerance : float
            Required error tolerance.
        criterion : {'default', 'RHS'}
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

    def _solve(self):

         if self.var.mesh.communicator.Nproc > 1:
             raise Exception("SciPy solvers cannot be used with multiple processors")

         self.var[:] = numerix.reshape(self._solve_(self.matrix, self.var.ravel(), numerix.array(self.RHSvector)), self.var.shape)
