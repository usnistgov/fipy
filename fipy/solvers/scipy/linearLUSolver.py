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

import os
import sys

from scipy.sparse.linalg import splu

from fipy.solvers.scipy.scipySolver import ScipySolver
from fipy.tools import numerix

DEBUG = False

class LinearLUSolver(ScipySolver):
    """
    
    The `LinearLUSolver` solves a linear system of equations using
    LU-factorisation. This method solves systems with a general
    non-symmetric coefficient matrix using partial pivoting.

    The `LinearLUSolver` is a wrapper class for the the PySparse_
    `superlu.factorize()` method.

    .. _PySparse: http://pysparse.sourceforge.net
    
    """
    
    def __init__(self, tolerance=1e-10, iterations=10, steps=None,
                       maxIterations=10, precon=None):
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
        
        super(LinearLUSolver, self).__init__(tolerance = tolerance, 
                                             iterations = iterations, 
                                             steps = steps)

    def _solve_(self, L, x, b):
        """
        :Parameters:
            - `L`: a `scipy.sparse.linalg.csr_matrix`.

        TODO: This method doesn't mirror its Pysparse counterpart in that `L` is
        a foreign sparse matrix instead of a FiPy sparse matrix.
        """
        diag = L.takeDiagonal()
        maxdiag = max(numerix.absolute(diag))

        L = L * (1 / maxdiag)
        b = b * (1 / maxdiag)

        import sys
        print >> sys.stderr, l.matrix.todense()

        LU = splu(L.matrix.asformat("csc"), diag_pivot_thresh=1.,
                                            drop_tol=0.,
                                            relax=1,
                                            panel_size=10,
                                            permc_spec=3)

        error0 = numerix.sqrt(numerix.sum((L * x - b)**2))

        for iteration in range(self.iterations):
            errorVector = L * x - b

            if (numerix.sqrt(numerix.sum(errorVector**2)) / error0)  <= self.tolerance:
                break

            xError = LU.solve(errorVector)
            x[:] = x - xError
            
        if os.environ.has_key('FIPY_VERBOSE_SOLVER'):
            from fipy.tools.debug import PRINT        
            PRINT('iterations: %d / %d' % (iteration+1, self.iterations))
            PRINT('residual:', numerix.sqrt(numerix.sum(errorVector**2)))
      
