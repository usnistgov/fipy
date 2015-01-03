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

            mesh = self.var.mesh
            localNonOverlappingCellIDs = mesh._localNonOverlappingCellIDs
            
            ## The following conditional is required because empty indexing is not altogether functional.
            ## This numpy.empty((0,))[[]] and this numpy.empty((0,))[...,[]] both work, but this
            ## numpy.empty((3, 0))[...,[]] is broken.
            if self.var.shape[-1] != 0:
                s = (Ellipsis, localNonOverlappingCellIDs)
            else:
                s = (localNonOverlappingCellIDs,)
                
            comm = self.var.mesh.communicator.petsc4py_comm
            nonOverlappingVector = PETSc.Vec().createWithArray(self.var[s].ravel(), comm=comm)
            
            from fipy.variables.coupledCellVariable import _CoupledCellVariable

            if isinstance(self.RHSvector, _CoupledCellVariable):
                RHSvector = self.RHSvector[localNonOverlappingCellIDs]
            else:
                RHSvector = numerix.reshape(numerix.array(self.RHSvector), self.var.shape)[s].ravel()
                
                    
            nonOverlappingRHSvector = PETSc.Vec().createWithArray(RHSvector)

            del RHSvector
            
            # PETSc requires that ghosts be at the end
            corporeals = numerix.in1d(mesh._globalOverlappingCellIDs, mesh._globalNonOverlappingCellIDs)
            ghosts = mesh._localOverlappingCellIDs[~corporeals]
            ids = numerix.concatenate([mesh._localOverlappingCellIDs[corporeals], ghosts])
            overlappingVector = PETSc.Vec().createGhostWithArray(ghosts=ghosts.astype('int32'), array=numerix.asarray(self.var[..., ids]), comm=PETSc.COMM_SELF)

            self.globalVectors = (globalMatrix, nonOverlappingVector, nonOverlappingRHSvector, overlappingVector)

        return self.globalVectors

    def _deleteGlobalMatrixAndVectors(self):
        self.matrix.flush()
        del self.globalVectors
        
    @property
    def global2local(self):
        if not hasattr(self, "_global2local"):
            globalMatrix, nonOverlappingVector, nonOverlappingRHSvector, overlappingVector = self._globalMatrixAndVectors
            
            fr = PETSc.IS().createGeneral(globalMatrix._globalOverlappingColIDs.astype('int32'), PETSc.COMM_SELF)
            to = PETSc.IS().createGeneral(globalMatrix._localOverlappingColIDs.astype('int32'), PETSc.COMM_SELF)
            self._global2local = PETSc.Scatter().create(nonOverlappingVector, fr, overlappingVector, to)
            fr.destroy()
            to.destroy()
        
        return self._global2local

    def _solve(self):
        from fipy.terms import SolutionVariableNumberError
        
        globalMatrix, nonOverlappingVector, nonOverlappingRHSvector, overlappingVector = self._globalMatrixAndVectors

        if ((self.matrix == 0)
            or (self.matrix.matrix.sizes[0][1] != self.matrix.matrix.sizes[1][1])
            or (self.matrix.matrix.sizes[0][1] != nonOverlappingVector.size)):

            raise SolutionVariableNumberError

        self._solve_(globalMatrix.matrix, 
                     nonOverlappingVector, 
                     nonOverlappingRHSvector)

        self.global2local.scatter(nonOverlappingVector, overlappingVector)
        
        with overlappingVector.localForm() as lf:
            self.var.value = numerix.reshape(numerix.array(lf), self.var.shape)
        
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
            
            mesh = self.var.mesh
            
            corporeals = numerix.in1d(mesh._globalOverlappingCellIDs, mesh._globalNonOverlappingCellIDs)
            ghosts = mesh._localOverlappingCellIDs[~corporeals]
            ids = numerix.concatenate([mesh._localOverlappingCellIDs[corporeals], ghosts])
            overlappingResidual = PETSc.Vec().createGhostWithArray(ghosts=ghosts.astype('int32'), array=numerix.asarray(self.var[..., ids]), comm=PETSc.COMM_SELF)
            
            self.global2local.scatter(residual, overlappingResidual)

            return overlappingResidual

    def _calcResidualVectorNonOverlapping_(self):
        globalMatrix, nonOverlappingVector, nonOverlappingRHSvector, overlappingVector = self._globalMatrixAndVectors
        residual = globalMatrix * nonOverlappingVector - nonOverlappingRHSvector
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
