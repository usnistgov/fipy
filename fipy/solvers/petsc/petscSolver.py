#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "trilinosSolver.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #  Author: Maxsim Gibiansky <maxsim.gibiansky@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #  
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this software is not subject to copyright
 # protection and is in the public domain.  FiPy is an experimental
 # system.  NIST assumes no responsibility whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 # 
 # This software can be redistributed and/or modified freely
 # provided that any derivative works bear some notice that they are
 # derived from it, and any modified versions bear some notice that
 # they have been modified.
 # ========================================================================
 #  
 # ###################################################################
 ##

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
            raise NotImplementedError, "can't instantiate abstract base class"
        else:
            Solver.__init__(self, *args, **kwargs)

    @property
    def _globalMatrixAndVectors(self):
        if not hasattr(self, 'globalVectors'):
            globalMatrix = self.matrix

            var = numerix.asarray(self.matrix._ghostTake(self.var))
            comm = self.var.mesh.communicator.petsc4py_comm
            overlappingVector = PETSc.Vec().createGhostWithArray(ghosts=self.matrix._ghosts.astype('int32'), 
                                                                 array=var.ravel(), 
                                                                 comm=comm)


            from fipy.variables.coupledCellVariable import _CoupledCellVariable
            if isinstance(self.RHSvector, _CoupledCellVariable):
                RHSvector = self.RHSvector
            else:
                RHSvector = numerix.reshape(numerix.asarray(self.RHSvector), self.var.shape)
                
            RHSvector = numerix.asarray(self.matrix._ghostTake(RHSvector))
            overlappingRHSvector = PETSc.Vec().createGhostWithArray(ghosts=self.matrix._ghosts.astype('int32'), 
                                                                    array=RHSvector.ravel(), 
                                                                    comm=comm)

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

        overlappingVector.ghostUpdate()
        with overlappingVector.localForm() as lf:
            value = numerix.reshape(numerix.array(lf), self.var.shape)
            self.matrix._ghostPut(self.var, value)
        
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
            residual, globalMatrix = self._calcResidualVectorNonOverlapping_()
            
            return residual

    def _calcResidualVectorNonOverlapping_(self):
        globalMatrix, overlappingVector, overlappingRHSvector = self._globalMatrixAndVectors
        residual = globalMatrix * overlappingVector - overlappingRHSvector
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
