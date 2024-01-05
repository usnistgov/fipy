from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from PyTrilinos import Epetra
from PyTrilinos import EpetraExt

from fipy.solvers.solver import Solver
from fipy.tools import numerix

class TrilinosSolver(Solver):

    """
    .. attention:: This class is abstract. Always create one of its subclasses.

    """

    def _storeMatrix(self, var, matrix, RHSvector):
        self.var = var
        if hasattr(self, 'matrix'):
            self.matrix.matrix = matrix.matrix
        else:
            self.matrix = matrix
        self.RHSvector = RHSvector

    @property
    def _globalMatrixAndVectors(self):
        if not hasattr(self, '_globalVectors'):
            globalMatrix = self.matrix.asTrilinosMeshMatrix()

            mesh = self.var.mesh
            localNonOverlappingCellIDs = mesh._localNonOverlappingCellIDs

            ## The following conditional is required because empty indexing is not altogether functional.
            ## This numpy.empty((0,))[[]] and this numpy.empty((0,))[...,[]] both work, but this
            ## numpy.empty((3, 0))[...,[]] is broken.
            if self.var.shape[-1] != 0:
                s = (Ellipsis, localNonOverlappingCellIDs)
            else:
                s = (localNonOverlappingCellIDs,)

            nonOverlappingVector = Epetra.Vector(globalMatrix.domainMap,
                                                 self.var[s].ravel())
            from fipy.variables.coupledCellVariable import _CoupledCellVariable

            if isinstance(self.RHSvector, _CoupledCellVariable):
                RHSvector = self.RHSvector[localNonOverlappingCellIDs]
            else:
                RHSvector = numerix.reshape(numerix.array(self.RHSvector), self.var.shape)[s].ravel()


            nonOverlappingRHSvector = Epetra.Vector(globalMatrix.rangeMap,
                                                    RHSvector)

            del RHSvector

            overlappingVector = Epetra.Vector(globalMatrix.colMap, self.var)

            self._globalVectors = (globalMatrix, nonOverlappingVector, nonOverlappingRHSvector, overlappingVector)

        return self._globalVectors

    def _deleteGlobalMatrixAndVectors(self):
        self.matrix.flush()
        del self._globalVectors

    def _rhsNorm(self, L, x, b):
        return float(b.Norm2())

    def _matrixNorm(self, L, x, b):
        return L.NormInf()

    def _residualVectorAndNorm(self, L, x, b):
        # residualVector = L*x - b
        residualVector = Epetra.Vector(L.RangeMap())
        L.Multiply(False, x, residualVector)
        # If A is an Epetra.Vector with map M
        # and B is an Epetra.Vector with map M
        # and C = A - B
        # then C is an Epetra.Vector with *no map* !!!?!?!
        residualVector -= b

        return residualVector, float(residualVector.Norm2())

    @property
    def _Lxb(self):
        """Matrix, solution vector, and right-hand side vector

        Returns
        -------
        L : Epetra.CrsMatrix
            Sparse matrix
        x : Epetra.Vector
            Solution variable as non-ghosted vector
        b : Epetra.Vector
            Right-hand side as non-ghosted vector
        """
        L, x, b, _ = self._globalMatrixAndVectors

        if not (L.rangeMap.SameAs(L.domainMap)
                and L.rangeMap.SameAs(x.Map())):

            from fipy.terms import SolutionVariableNumberError

            raise SolutionVariableNumberError

        return (L.matrix, x, b)

    def _scatterGhosts(self, x):
        """Distribute ghost values (if any) across processes
        """
        globalMatrix, _, _, overlappingVector = self._globalMatrixAndVectors

        overlappingVector.Import(x,
                                 Epetra.Import(globalMatrix.colMap,
                                               globalMatrix.domainMap),
                                 Epetra.Insert)

        return numerix.asarray(overlappingVector)

    def _cleanup(self):
        self._deleteGlobalMatrixAndVectors()
        del self.var
        del self.RHSvector

    @property
    def _matrixClass(self):
        from fipy.solvers import _MeshMatrix
        # could be Trilinos or Pysparse, depending on configuration
        return _MeshMatrix

    def _calcResidualVector(self, residualFn=None):
        if residualFn is not None:
            return residualFn(self.var, self.matrix, self.RHSvector)
        else:
            residual, globalMatrix = self._calcResidualVectorNonOverlapping_()

            overlappingResidual = Epetra.Vector(globalMatrix.colMap)
            overlappingResidual.Import(residual,
                                       Epetra.Import(globalMatrix.colMap,
                                                     globalMatrix.domainMap),
                                       Epetra.Insert)

            return overlappingResidual

    def _calcResidualVectorNonOverlapping_(self):
        globalMatrix, nonOverlappingVector, nonOverlappingRHSvector, overlappingVector = self._globalMatrixAndVectors
        # If A is an Epetra.Vector with map M
        # and B is an Epetra.Vector with map M
        # and C = A - B
        # then C is an Epetra.Vector with *no map* !!!?!?!
        residual = globalMatrix * nonOverlappingVector
        residual -= nonOverlappingRHSvector
        return residual, globalMatrix

    def _calcResidual(self, residualFn=None):
        if residualFn is not None:
            return residualFn(self.var, self.matrix, self.RHSvector)
        else:
            comm = self.var.mesh.communicator
            residual, globalMatrix = self._calcResidualVectorNonOverlapping_()
            return comm.Norm2(residual)

    def _calcRHSNorm(self):
        return self.nonOverlappingRHSvector.Norm2()
