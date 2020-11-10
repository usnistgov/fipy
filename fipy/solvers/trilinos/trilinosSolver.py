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
    def __init__(self, *args, **kwargs):
        if self.__class__ is TrilinosSolver:
            raise NotImplementedError("can't instantiate abstract base class")
        else:
            Solver.__init__(self, *args, **kwargs)

    def _storeMatrix(self, var, matrix, RHSvector):
        self.var = var
        if hasattr(self, 'matrix'):
            self.matrix.matrix = matrix.matrix
        else:
            self.matrix = matrix
        self.RHSvector = RHSvector

    @property
    def _globalMatrixAndVectors(self):
        if not hasattr(self, 'globalVectors'):
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

            self.globalVectors = (globalMatrix, nonOverlappingVector, nonOverlappingRHSvector, overlappingVector)

        return self.globalVectors

    def _deleteGlobalMatrixAndVectors(self):
        self.matrix.flush()
        del self.globalVectors

    def _solve(self):
        from fipy.terms import SolutionVariableNumberError

        globalMatrix, nonOverlappingVector, nonOverlappingRHSvector, overlappingVector = self._globalMatrixAndVectors

        if not (globalMatrix.rangeMap.SameAs(globalMatrix.domainMap)
                and globalMatrix.rangeMap.SameAs(nonOverlappingVector.Map())):

            raise SolutionVariableNumberError

        self._solve_(globalMatrix.matrix,
                     nonOverlappingVector,
                     nonOverlappingRHSvector)

        overlappingVector.Import(nonOverlappingVector,
                                 Epetra.Import(globalMatrix.colMap,
                                               globalMatrix.domainMap),
                                 Epetra.Insert)

        self.var.value = numerix.reshape(numerix.array(overlappingVector), self.var.shape)

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
