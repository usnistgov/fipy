#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "trilinosAztecOOSolver.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
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

import os

from petsc4py import PETSc

from fipy.solvers.petsc.petscSolver import PETScSolver

__all__ = ["LinearLUSolver"]

class LinearLUSolver(PETScSolver):

    """
    The `LinearLUSolver` is an interface to the LU preconditioner in PETSc.
    A direct solve is performed.

    """
      
    def __init__(self, tolerance=1e-10, iterations=1000, precon="lu"):
        """
        :Parameters:
          - `tolerance`: The required error tolerance.
          - `iterations`: The maximum number of iterative steps to perform.
          - `precon`: *Ignored*. 

        """
        PETScSolver.__init__(self, tolerance=tolerance,
                             iterations=iterations, precon="lu")

    def _solve_(self, L, x, b):
        ksp = PETSc.KSP()
        ksp.create(PETSc.COMM_WORLD)
        ksp.setType("preonly")
        ksp.getPC().setType(self.preconditioner)
        # TODO: SuperLU invoked with PCFactorSetMatSolverPackage(pc, MATSOLVERSUPERLU)
        #       see: http://www.mcs.anl.gov/petsc/petsc-dev/src/ksp/ksp/examples/tutorials/ex52.c.html
        # PETSc.PC().setFactorSolverPackage("superlu")
        
        L.assemblyBegin()
        L.assemblyEnd()
        ksp.setOperators(L)
        ksp.setFromOptions()
        
        for iteration in range(self.iterations):
            errorVector = L * x - b
            tol = errorVector.norm()
            
            if iteration == 0:
                tol0 = tol
                
            if (tol / tol0) <= self.tolerance:
                break
                
            xError = x.copy()

            ksp.solve(errorVector, xError)
            x -= xError
            
        if 'FIPY_VERBOSE_SOLVER' in os.environ:
            from fipy.tools.debug import PRINT        
            PRINT('iterations: %d / %d' % (iteration+1, self.iterations))
            PRINT('residual:', errorVector.narm(1))

