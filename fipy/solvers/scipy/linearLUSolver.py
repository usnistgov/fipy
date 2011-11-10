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

from scipy.sparse.linalg import splu

from fipy.solvers.scipy.scipySolver import _ScipySolver
from fipy.tools import numerix

__all__ = ["LinearLUSolver"]

class LinearLUSolver(_ScipySolver):
    """
    The `LinearLUSolver` solves a linear system of equations using
    LU-factorisation.  The `LinearLUSolver` is a wrapper class for the
    the Scipy `scipy.sparse.linalg.splu` moduleq.
    """
    
    def _solve_(self, L, x, b):
        diag = L.takeDiagonal()
        maxdiag = max(numerix.absolute(diag))

        L = L * (1 / maxdiag)
        b = b * (1 / maxdiag)

        LU = splu(L.matrix.asformat("csc"), diag_pivot_thresh=1.,
                                            drop_tol=0.,
                                            relax=1,
                                            panel_size=10,
                                            permc_spec=3)

        error0 = numerix.sqrt(numerix.sum((L * x - b)**2))

        for iteration in range(min(self.iterations, 10)):
            errorVector = L * x - b

            if (numerix.sqrt(numerix.sum(errorVector**2)) / error0)  <= self.tolerance:
                break

            xError = LU.solve(errorVector)
            x[:] = x - xError
            
        if 'FIPY_VERBOSE_SOLVER' in os.environ:
            from fipy.tools.debug import PRINT        
            PRINT('iterations: %d / %d' % (iteration+1, self.iterations))
            PRINT('residual:', numerix.sqrt(numerix.sum(errorVector**2)))

        return x
