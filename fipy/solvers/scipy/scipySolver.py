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

    @property
    def _matrixClass(self):
        return _ScipyMeshMatrix

    def _solve(self):

         if self.var.mesh.communicator.Nproc > 1:
             raise Exception("SciPy solvers cannot be used with multiple processors")

         self.var[:] = numerix.reshape(self._solve_(self.matrix, self.var.ravel(), numerix.array(self.RHSvector)), self.var.shape)
