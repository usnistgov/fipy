__docformat__ = 'restructuredtext'

from petsc4py import PETSc

from fipy.solvers.solver import Solver
from fipy.tools import numerix
from fipy.matrices.petscMatrix import _PETScMeshMatrix

class PETScSolver(Solver):

    """
    .. attention:: This class is abstract. Always create one of its subclasses.

    """

    @property
    def globalVectors(self):
        if not hasattr(self, '_globalVectors'):
            overlappingVector = self.matrix._fipy2petscGhost(var=self.var)

            from fipy.variables.coupledCellVariable import _CoupledCellVariable
            if isinstance(self.RHSvector, _CoupledCellVariable):
                RHSvector = self.RHSvector
            else:
                RHSvector = numerix.reshape(numerix.asarray(self.RHSvector), self.var.shape)
                
            overlappingRHSvector = self.matrix._fipy2petscGhost(var=RHSvector)

            self._globalVectors = (overlappingVector, overlappingRHSvector)

        return self._globalVectors

    def _deleteGlobalVectors(self):
        if hasattr(self, "_globalVectors"):
            overlappingVector, overlappingRHSvector = self.globalVectors
            overlappingVector.destroy()
            overlappingRHSvector.destroy()
        del self._globalVectors
        
    def _rhsNorm(self, L, x, b):
        return b.norm(PETSc.NormType.NORM_2)

    def _matrixNorm(self, L, x, b):
        L.assemble()
        return L.norm(PETSc.NormType.NORM_INFINITY)

    def _residualVectorAndNorm(self, L, x, b):
        residualVector = L * x - b

        return residualVector, residualVector.norm(PETSc.NormType.NORM_2)

    @property
    def _Lxb(self):
        """Matrix, solution vector, and right-hand side vector

        Returns
        -------
        L : PETSc.Mat
            Sparse matrix
        x : PETSc.Vec
            Solution variable as ghosted vector
        b : PETSc.Vec
            Right-hand side as ghosted vector
        """
        x, b = self.globalVectors
        L = self.matrix.matrix

        if ((self.matrix == 0)
            or (L.sizes[0][1] != L.sizes[1][1])
            or (L.sizes[0][1] != x.size)):

            from fipy.terms import SolutionVariableNumberError

            raise SolutionVariableNumberError

        return (L, x, b)

    def _scatterGhosts(self, x):
        """Distribute ghost values (if any) across processes
        """
        return self.matrix._petsc2fipyGhost(vec=x)

    def _cleanup(self):
        self._deleteGlobalVectors()
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
        overlappingVector, overlappingRHSvector = self.globalVectors
        Lx = self.matrix * overlappingVector
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
        if hasattr(self, "_globalVectors"):
            overlappingVector, overlappingRHSvector = self.globalVectors
            del self.matrix
            overlappingVector.destroy()
            overlappingRHSvector.destroy()
