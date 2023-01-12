__docformat__ = 'restructuredtext'

from petsc4py import PETSc

from fipy.solvers.solver import Solver
from fipy.tools import numerix
from fipy.matrices.petscMatrix import _PETScMeshMatrix

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
        del self.matrix
        if hasattr(self, "globalVectors"):
            globalMatrix, overlappingVector, overlappingRHSvector = self._globalMatrixAndVectors
            overlappingVector.destroy()
            overlappingRHSvector.destroy()
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
        return _PETScMeshMatrix

    def _calcResidualVector(self, residualFn=None):
        if residualFn is not None:
            return residualFn(self.var, self.matrix, self.RHSvector)
        else:
            residual = self._calcResidualVector_()
            
            residual.ghostUpdate()
            with residual.localForm() as lf:
                arr = numerix.array(lf)
            residual.destroy()
            return arr

    def _calcResidualVector_(self):
        globalMatrix, overlappingVector, overlappingRHSvector = self._globalMatrixAndVectors
        Lx = globalMatrix * overlappingVector
        residual = Lx - overlappingRHSvector
        Lx.destroy()
        return residual

    def _calcResidual(self, residualFn=None):
        if residualFn is not None:
            return residualFn(self.var, self.matrix, self.RHSvector)
        else:
            comm = self.var.mesh.communicator
            residual = self._calcResidualVector_()
            norm = comm.Norm2(residual)
            residual.destroy()
            return norm
        
    def _calcRHSNorm(self):
        return self.nonOverlappingRHSvector.Norm2()

    def __del__(self):
        if hasattr(self, "globalVectors"):
            globalMatrix, overlappingVector, overlappingRHSvector = self._globalMatrixAndVectors
            del globalMatrix
            overlappingVector.destroy()
            overlappingRHSvector.destroy()
