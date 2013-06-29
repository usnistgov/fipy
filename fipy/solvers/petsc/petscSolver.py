__docformat__ = 'restructuredtext'

from fipy.solvers.solver import Solver
from fipy.tools import numerix

class PETScSolver(Solver):

    """
    .. attention:: This class is abstract. Always create one of its subclasses.

    """
    def __init__(self, *args, **kwargs):
        if self.__class__ is PETScSolver:
            raise NotImplementedError, "can't instantiate abstract base class"
        else:
            Solver.__init__(self, *args, **kwargs)

    def _solve(self):
        array = self.var.numericValue.ravel()

        from fipy.terms import SolutionVariableNumberError
        
        if ((self.matrix == 0)
            or (self.matrix._shape[0] != self.matrix._shape[1])
            or (self.matrix._shape[0] != len(array))):

            raise SolutionVariableNumberError
        
        self._solve_(self.matrix, array, self.RHSvector)
        factor = self.var.unit.factor
        if factor != 1:
            array /= self.var.unit.factor

        self.var[:] = array.reshape(self.var.shape)
            
    @property
    def _matrixClass(self):
        from fipy.solvers import _MeshMatrix
        return _MeshMatrix
