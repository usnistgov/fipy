__docformat__ = 'restructuredtext'

from petsc4py import PETSc

from fipy.solvers.solver import Solver
from fipy.tools import numerix

class PETScSolver(Solver):

    """
    .. attention:: This class is abstract. Always create one of its subclasses.

    """
    def __init__(self, *args, **kwargs):
        if self.__class__ is PETScSolver:
            raise NotImplementedError("can't instantiate abstract base class")
        else:
            Solver.__init__(self, *args, **kwargs)

    @property
    def _globalMatrixAndVectors(self):
        if not hasattr(self, 'globalVectors'):
            globalMatrix = self.matrix

            overlappingVector = self.matrix._fipy2petscGhost(var=self.var)

            from fipy.variables.coupledCellVariable import _CoupledCellVariable
            if isinstance(self.RHSvector, _CoupledCellVariable):
                RHSvector = self.RHSvector
            else:
                RHSvector = numerix.reshape(numerix.asarray(self.RHSvector), self.var.shape)
                
            overlappingRHSvector = self.matrix._fipy2petscGhost(var=RHSvector)

            self.globalVectors = (globalMatrix, overlappingVector, overlappingRHSvector)

        return self.globalVectors

    def _deleteGlobalMatrixAndVectors(self):
        self.matrix.flush()
        del self.globalVectors
        
    def _solve(self):
        from fipy.terms import SolutionVariableNumberError
        
        globalMatrix, overlappingVector, overlappingRHSvector = self._globalMatrixAndVectors

        if ((self.matrix == 0)
            or (self.matrix.matrix.sizes[0][1] != self.matrix.matrix.sizes[1][1])
            or (self.matrix.matrix.sizes[0][1] != overlappingVector.size)):

            raise SolutionVariableNumberError

        self._solve_(globalMatrix.matrix, 
                     overlappingVector, 
                     overlappingRHSvector)

        value = self.matrix._petsc2fipyGhost(vec=overlappingVector)
        self.var.value = numerix.reshape(value, self.var.shape)
        
        self._deleteGlobalMatrixAndVectors()
        del self.var
        del self.RHSvector
            
    @property
    def _matrixClass(self):
        from fipy.solvers import _MeshMatrix
        return _MeshMatrix

    def _calcResidualVector(self, residualFn=None):
        if residualFn is not None:
            return residualFn(self.var, self.matrix, self.RHSvector)
        else:
            residual, globalMatrix = self._calcResidualVector_()
            
            residual.ghostUpdate()
            with residual.localForm() as lf:
                residual = numerix.array(lf)
            return residual

    def _calcResidualVector_(self):
        globalMatrix, overlappingVector, overlappingRHSvector = self._globalMatrixAndVectors
        residual = globalMatrix * overlappingVector - overlappingRHSvector
        return residual, globalMatrix

    def _calcResidual(self, residualFn=None):
        if residualFn is not None:
            return residualFn(self.var, self.matrix, self.RHSvector)
        else:
            comm = self.var.mesh.communicator
            residual, globalMatrix = self._calcResidualVector_()
            return comm.Norm2(residual)
        
    def _calcRHSNorm(self):
        return self.nonOverlappingRHSvector.Norm2()
