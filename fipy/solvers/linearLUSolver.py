#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "linearLUSolver.py"
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

import sys

import precon
import itsolvers
import superlu

from fipy.solvers.solver import _Solver
from fipy.tools.sparseMatrix import _SparseMatrix
import Numeric

class LinearLUSolver(_Solver):
    """
    
    The `LinearLUSolver` solves a linear system of equations using
    LU-factorisation. This method solves systems with a general
    non-symmetric coefficient matrix using partial pivoting.

    The `LinearLUSolver` is a wrapper class for the the PySparse_
    `superlu.factorize` method. Usage:

    ::

        solver = LinearLUSolver(tolerance = 1e-10, steps = 1000)

    .. _PySparse: http://pysparse.sourceforge.net
    
    """

    def __init__(self, tolerance = 1e-10, steps = 10):
        """
        Creates a `LinearLUSolver`.

        :Parameters:
          - `tolerance` : The required error tolerance.
          - `steps` : The number of LU decompositions to perform.
            For large systems a number of steps is generally required.

        """
        
	_Solver.__init__(self, tolerance = tolerance, steps = steps)

    def _solve(self, L, x, b):
        LU = superlu.factorize(L._getMatrix().to_csr(), diag_pivot_thresh = 1., relax = 1)
        LU.solve(b, x)

        
        tol = self.tolerance + 1.
        step = 0

        while tol > self.tolerance and step < self.steps:

            errorVector = L * x - b
            xError = Numeric.zeros(len(b),'d')
            LU = superlu.factorize(L._getMatrix().to_csr())
            LU.solve(errorVector, xError)
            x[:] = x - xError

            arg = Numeric.argmax(Numeric.absolute(xError))
            tol = Numeric.absolute(xError)[arg]
            step += 1

