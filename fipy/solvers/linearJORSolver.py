#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "linearJORSolver.py"
 #                                    created: 11/14/03 {3:56:49 PM} 
 #                                last update: 10/26/04 {11:37:29 AM} 
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
 #  Description: 
 # 
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2003-11-14 JEG 1.0 original
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

from fipy.solvers.solver import _Solver
from fipy.tools.sparseMatrix import _SparseMatrix

class LinearJORSolver(_Solver):
    """
    
    The `LinearJORSolver` solves a linear system of equations using
    Jacobi over-relaxation. This method solves systems with a general
    non-symmetric coefficient matrix.

    `superlu.factorize` method. Usage:

    ::

        solver = LinearJORSolver(tolerance = 1e-10, steps = 1000, relaxation = 1.0)

    """
    def __init__(self, tolerance = 1e-10, steps = 1000, relaxation = 1.0):
        """
        The `Solver` class should not be invoked directly.

        :Parameters:
          - `tolerance`: The required error tolerance.
          - `steps`: The maximum number of iterative steps to perform.
          - `relaxation`: The relaxation.
          
        """
	self.tolerance = tolerance
	self.steps = steps
        self.relaxation = relaxation
        
    def _solve(self, L, x, b):

        d = L.takeDiagonal()
        D = _SparseMatrix(len(d))
        D.putDiagonal(d)

        LU = L - D
        step = 0
        tol = 1e+10
        xold = x.copy()

        while step < self.steps and tol > self.tolerance:
            residual = L * x - b

            xold[:] = x
            x[:] = (-(LU) * x + b) / d

            x[:] = xold + self.relaxation * (x - xold)  

            tol = max(abs(residual))

            step += 1
            
            print step,tol
            
