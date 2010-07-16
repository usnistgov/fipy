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

import sys

from fipy.solvers.trilinos.trilinosSolver import TrilinosSolver
from fipy.matrices.trilinosMatrix import _trilinosToNumpyVector
from fipy.matrices.trilinosMatrix import _numpyToTrilinosVector

from fipy.tools import numerix

from PyTrilinos import Epetra
from PyTrilinos import EpetraExt
from PyTrilinos import Amesos

class LinearLUSolver(TrilinosSolver):

    """
    The `LinearLUSolver` is an interface to the Amesos KLU solver in Trilinos.

    """

    def __init__(self, tolerance=1e-10, iterations=10, steps=None, precon=None, maxIterations=10):
        """
        :Parameters:
          - `tolerance`: The required error tolerance.
          - `iterations`: The maximum number of iterative steps to perform.
          - `steps`: A deprecated name for `iterations`.

        """

        iterations = min(iterations, maxIterations)
        
        TrilinosSolver.__init__(self, tolerance=tolerance, 
                                iterations=iterations, steps=steps, precon=None)

        if precon is not None:
            import warnings
            warnings.warn("Trilinos KLU solver does not accept preconditioners.",
                           UserWarning, stacklevel=2)
        self.Factory = Amesos.Factory()

       
    def _solve_(self, L, x, b):
         
        for iteration in range(self.iterations):

             # errorVector = L*x - b
             errorVector = Epetra.Vector(L.RowMap())
             L.Multiply(False, x, errorVector)
             errorVector = errorVector - b

             tol = max(numerix.absolute(_trilinosToNumpyVector(errorVector)))

             if iteration == 0:
                 tol0 = tol

             if (tol / tol0) <= self.tolerance: 
                 break

             xError = _numpyToTrilinosVector(numerix.zeros(errorVector.GlobalLength(), 'd'), L.RowMap())
              
             Problem = Epetra.LinearProblem(L, xError, errorVector)
             Solver = self.Factory.Create("Klu", Problem)
             Solver.Solve()

             x[:] = x - xError
