#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "linearLUSolver.py"
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

import sys

from pysparse import precon
from pysparse import itsolvers
from pysparse import superlu

from fipy.solvers.pysparse.pysparseSolver import PysparseSolver
from fipy.tools import numerix

class LinearLUSolver(PysparseSolver):
    """
    
    The `LinearLUSolver` solves a linear system of equations using
    LU-factorisation. This method solves systems with a general
    non-symmetric coefficient matrix using partial pivoting.

    The `LinearLUSolver` is a wrapper class for the the PySparse_
    `superlu.factorize()` method.

    .. _PySparse: http://pysparse.sourceforge.net
    
    """
    
    def __init__(self, tolerance=1e-10, iterations=10, steps=None, precon=None, maxIterations=10):
        """
        Creates a `LinearLUSolver`.

        :Parameters:
          - `tolerance`: The required error tolerance.
          - `iterations`: The number of LU decompositions to perform.
          - `steps`: A deprecated name for `iterations`.
            For large systems a number of iterations is generally required.
          - `precon`: not used but maintains a common interface.

        """

        iterations = min(iterations, maxIterations)
        
        PysparseSolver.__init__(self, tolerance = tolerance, iterations=iterations, steps = steps, precon = precon)

    def _solve_(self, L, x, b):
        diag = L.takeDiagonal()
        maxdiag = max(numerix.absolute(diag))

        L = L * (1 / maxdiag)
        b = b * (1 / maxdiag)

        LU = superlu.factorize(L._getMatrix().to_csr())

        error0 = numerix.sqrt(numerix.sum((L * x - b)**2))

        for iteration in range(self.iterations):
            errorVector = L * x - b

            if (numerix.sqrt(numerix.sum(errorVector**2)) / error0)  <= self.tolerance:
                break

            xError = numerix.zeros(len(b),'d')
            LU.solve(errorVector, xError)
            x[:] = x - xError
