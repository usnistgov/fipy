#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "linearScipyLUSolver.py"
 #                                    created: 11/14/03 {3:56:49 PM} 
 #                                last update: 9/3/04 {10:43:22 PM} 
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

import Numeric
import scipy.linalg

from fipy.solvers.solver import Solver

class LinearScipyLUSolver(Solver):
    """
    The `LinearScipyLUSolver` solves a system of equations using
    LU-factorisation. This method solves systems of general
    non-symmetric matrices with partial pivoting.

    The `LinearScipyLUSolver` is a wrapper class for the the Scipy_
    `scipy.linalg.lu_solve` method. Usage:

    ::

        solver = LinearScipyLUSolver(tolerance = 1e-10, steps = 1000)

    .. _Scipy: http://www.scipy.org

    """
    
    def __init__(self, tolerance = 1e-10, steps = 10):
        """
        Creates a `LinearScipyLUSolver`.

        :Parameters:
          - `tolerance` : The required error tolerance.
          - `steps` : The number of LU decompositions to perform.
            For large systems a number of steps is generally required.

        """
	Solver.__init__(self, tolerance = tolerance, steps = steps)

    def _solve(self, L, x, b):

##        x[:] = scipy.linalg.solve(Numeric.array(L), b)
        LU = scipy.linalg.lu_factor(Numeric.array(L))
        x[:] = scipy.linalg.lu_solve(LU, b)
        tol = self.tolerance + 1.
        step = 0

        while tol > self.tolerance and step < self.steps:

            errorVector = L * x - b
            LU = scipy.linalg.lu_factor(Numeric.array(L))
            xError = scipy.linalg.lu_solve(LU, errorVector)
            x[:] = x - xError

            arg = Numeric.argmax(Numeric.absolute(xError))
            tol = Numeric.absolute(xError)[arg]
            step += 1
            
            print tol
