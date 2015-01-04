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
    def _bodies(self):
        if not hasattr(self, "_bodies_"):
            mesh = self.var.mesh
            self._bodies_ = numerix.in1d(mesh._globalOverlappingCellIDs, 
                                         mesh._globalNonOverlappingCellIDs)
        return self._bodies_
        
    @property
    def _ghosts(self):
        if not hasattr(self, "_ghosts_"):
            mesh = self.var.mesh
            self._ghosts_ = mesh._globalOverlappingCellIDs[~self._bodies]
            
        return self._ghosts_

    @property
    def _ghostSlice(self):
        if not hasattr(self, "_ghostSlice_"):
            mesh = self.var.mesh
            
            # PETSc requires that ghosts be at the end, FiPy doesn't care
            ids = numerix.concatenate([mesh._localOverlappingCellIDs[self._bodies], 
                                       mesh._localOverlappingCellIDs[~self._bodies]])
            
            ## The following conditional is required because empty indexing is not altogether functional.
            ## This numpy.empty((0,))[[]] and this numpy.empty((0,))[...,[]] both work, but this
            ## numpy.empty((3, 0))[...,[]] is broken.
            if self.var.shape[-1] != 0:
                self._ghostSlice_ = (Ellipsis, ids)
            else:
                self._ghostSlice_ = (ids,)


            # g                 g
            # 0 1 2 3 4 5 6 7 8 9  FiPy
            
            #                 g g
            # 1 2 3 4 5 6 7 8 0 9  PETSc
            
            # 8 0 1 2 3 4 5 6 7 9
        return self._ghostSlice_

    @property
    def _globalMatrixAndVectors(self):
        if not hasattr(self, 'globalVectors'):
            globalMatrix = self.matrix

            overlappingVector = PETSc.Vec().createGhostWithArray(ghosts=self._ghosts.astype('int32'), 
                                                                 array=numerix.asarray(self.var[self._ghostSlice]).ravel(), 
                                                                 comm=PETSc.COMM_WORLD)

            overlappingRHSvector = PETSc.Vec().createGhostWithArray(ghosts=self._ghosts.astype('int32'), 
                                                                    array=numerix.asarray(self.RHSvector[self._ghostSlice]).ravel(), 
                                                                    comm=PETSc.COMM_WORLD)

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
            self.var.value[self._ghostSlice] = numerix.reshape(numerix.array(lf), self.var.shape)
        
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
            
            residual.ghostUpdate()
            with residual.localForm() as lf:
                overlappingResidual = self.var.copy()
                overlappingResidual[self._ghostSlice] = numerix.reshape(numerix.array(lf), self.var.shape)

            return overlappingResidual

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
