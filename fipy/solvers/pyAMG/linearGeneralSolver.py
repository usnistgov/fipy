#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "linearGeneralSolver.py"
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
 
from fipy.solvers.scipy.scipySolver import _ScipySolver
from pyamg import solve
import os
from fipy.tools import numerix

__all__ = ["LinearGeneralSolver"]

class LinearGeneralSolver(_ScipySolver):
    """
    The `LinearGeneralSolver` is an interface to the generic pyAMG,
    which solves the arbitrary system Ax=b with the best out-of-the box
    choice for a solver. See `pyAMG.solve` for details.
    """

    def _solve_(self, L, x, b):
        if 'FIPY_VERBOSE_SOLVER' in os.environ:
            verbosity = True
        else:
            verbosity = False

        return solve(L.matrix, b, verb=verbosity, tol=self.tolerance)

