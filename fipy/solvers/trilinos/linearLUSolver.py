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
 # of Standards and Technology, an agency of the Federal Government.
 # Pursuant to title 17 section 105 of the United States Code,
 # United States Code this software is not subject to copyright
 # protection, and this software is considered to be in the public domain.
 # FiPy is an experimental system.
 # NIST assumes no responsibility whatsoever for its use by whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 #
 # To the extent that NIST may hold copyright in countries other than the
 # United States, you are hereby granted the non-exclusive irrevocable and
 # unconditional right to print, publish, prepare derivative works and
 # distribute this software, in any medium, or authorize others to do so on
 # your behalf, on a royalty-free basis throughout the world.
 #
 # You may improve, modify, and create derivative works of the software or
 # any portion of the software, and you may copy and distribute such
 # modifications or works.  Modified works should carry a notice stating
 # that you changed the software and should note the date and nature of any
 # such change.  Please explicitly acknowledge the National Institute of
 # Standards and Technology as the original source.
 #
 # This software can be redistributed and/or modified freely provided that
 # any derivative works bear some notice that they are derived from it, and
 # any modified versions bear some notice that they have been modified.
 # ========================================================================
 #
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

import os

from PyTrilinos import Epetra
from PyTrilinos import Amesos

from fipy.solvers.trilinos.trilinosSolver import TrilinosSolver

__all__ = ["LinearLUSolver"]

class LinearLUSolver(TrilinosSolver):

    """
    The `LinearLUSolver` is an interface to the Amesos KLU solver in Trilinos.

    """

    def __init__(self, tolerance=1e-10, iterations=10, precon=None, maxIterations=10):
        """
        :Parameters:
          - `tolerance`: The required error tolerance.
          - `iterations`: The maximum number of iterative steps to perform.

        """

        iterations = min(iterations, maxIterations)

        TrilinosSolver.__init__(self, tolerance=tolerance,
                                iterations=iterations, precon=None)

        if precon is not None:
            import warnings
            warnings.warn("Trilinos KLU solver does not accept preconditioners.",
                           UserWarning, stacklevel=2)
        self.Factory = Amesos.Factory()


    def _solve_(self, L, x, b):

        for iteration in range(self.iterations):
             # errorVector = L*x - b
             errorVector = Epetra.Vector(L.RangeMap())
             L.Multiply(False, x, errorVector)
             # If A is an Epetra.Vector with map M
             # and B is an Epetra.Vector with map M
             # and C = A - B
             # then C is an Epetra.Vector with *no map* !!!?!?!
             errorVector -= b

             tol = errorVector.Norm1()

             if iteration == 0:
                 tol0 = tol

             if (tol / tol0) <= self.tolerance:
                 break

             xError = Epetra.Vector(L.RowMap())

             Problem = Epetra.LinearProblem(L, xError, errorVector)
             Solver = self.Factory.Create("Klu", Problem)
             Solver.Solve()

             x[:] = x - xError

        if 'FIPY_VERBOSE_SOLVER' in os.environ:
            from fipy.tools.debug import PRINT
            PRINT('iterations: %d / %d' % (iteration + 1, self.iterations))
            PRINT('residual:', errorVector.Norm2())
