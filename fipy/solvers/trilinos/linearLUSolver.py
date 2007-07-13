#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "linearLUSolver.py"
 #                                    created: 06/07/07 
 #                                last update: 06/11/07 
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
 #  Description: 
 # 
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2007-06-12 MLG 1.0 original
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

import sys

from fipy.solvers.trilinos.trilinosSolver import TrilinosSolver

from PyTrilinos import Epetra
from PyTrilinos import Amesos

class LinearLUSolver(TrilinosSolver):

    """
    An interface to the Amesos KLU solver in Trilinos.

    """
    def __init__(self, tolerance=0, iterations=0, steps=None, precon=None):
        """
        :Parameters:
        - `tolerance`: The required error tolerance.
        - `iterations`: The maximum number of iterative steps to perform.
        - `steps`: A deprecated name for `iterations`.
        """
        TrilinosSolver.__init__(self, tolerance=tolerance, 
                                iterations=iterations, steps=steps, precon=None)

        if  tolerance !=0 or iterations != 0:
            import warnings
            warnings.warn("Trilinos KLU solver currently ignores tolerance and iteration specifications and runs a single iteration.", UserWarning, stacklevel=2)

        if precon is not None:
            import warnings
            warnings.warn("Trilinos KLU solver does not accept preconditioners.",
                           UserWarning, stacklevel=2)
        self.Factory = Amesos.Factory()

       
    def _applyTrilinosSolver(self, A, LHS, RHS):
        Problem = Epetra.LinearProblem(A, LHS, RHS)
        Solver = self.Factory.Create("Klu", Problem)
        Solver.Solve()
