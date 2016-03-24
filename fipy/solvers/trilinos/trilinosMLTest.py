#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "trilinosMLTest.py"
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

from PyTrilinos import Epetra
from PyTrilinos import EpetraExt
from PyTrilinos import Amesos
from PyTrilinos import AztecOO
from PyTrilinos import ML
from PyTrilinos import IFPACK

from fipy.solvers.trilinos.trilinosSolver import TrilinosSolver
from fipy.tools import numerix

__all__ = ["TrilinosMLTest"]

class TrilinosMLTest(TrilinosSolver):

    """
    This solver class does not actually solve the system, but outputs
    information about what ML preconditioner settings will work best.
    """

    def __init__(self, tolerance=1e-10, iterations=5, MLOptions={}, testUnsupported = False):
        """
        :Parameters:
          - `tolerance`: The required error tolerance.
          - `iterations`: The maximum number of iterations to perform per test.
          - `MLOptions`: Options to pass to ML. A dictionary of {option:value} pairs. This will be passed to ML.SetParameterList.
          - `testUnsupported`: test smoothers that are not currently implemented in preconditioner objects.

        For detailed information on the possible parameters for ML, see
        http://trilinos.sandia.gov/packages/ml/documentation.html

        Currently, passing options to Aztec through ML is not supported.
         """

        TrilinosSolver.__init__(self, tolerance=tolerance,
                                iterations=iterations)

        self.MLOptions = MLOptions
        if "output" not in self.MLOptions:
            self.MLOptions["output"] = 0

        if "test: max iters" not in self.MLOptions:
            self.MLOptions["test: max iters"] = iterations

        if "test: tolerance" not in self.MLOptions:
            self.MLOptions["test: tolerance"] = tolerance


        unsupportedSmoothers = ["Jacobi", "Gauss-Seidel", "block Gauss-Seidel", "ParaSails", "IFPACK", "ML"]

        if not testUnsupported:
            for smoother in unsupportedSmoothers:
                if ("test: " + smoother) not in self.MLOptions:
                    self.MLOptions["test: " + smoother] = False



    def _applyTrilinosSolver(self, A, LHS, RHS):

        Prec = ML.MultiLevelPreconditioner(A, False)

        Prec.SetParameterList(self.MLOptions)
        Prec.ComputePreconditioner()

        Prec.TestSmoothers()
        input("Results of preconditioner tests shown above. Currently, the first tests in the 'Gauss-Seidel (sym)','Aztec preconditioner', and 'Aztec as solver' sections indicate the expected performance of the MultilevelSGSPreconditioner, MultilevelDDPreconditioner, and MultilevelSolverSmootherPreconditioner classes, respectively.\n\nPress enter to quit.")
        import sys
        sys.exit(0)
