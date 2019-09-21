from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from PyTrilinos import Epetra
from PyTrilinos import EpetraExt
from PyTrilinos import Amesos
from PyTrilinos import AztecOO
from PyTrilinos import ML
from PyTrilinos import IFPACK

from fipy import input
from fipy.solvers.trilinos.trilinosSolver import TrilinosSolver
from fipy.tools import numerix

__all__ = ["TrilinosMLTest"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class TrilinosMLTest(TrilinosSolver):

    """
    This solver class does not actually solve the system, but outputs
    information about what ML preconditioner settings will work best.
    """

    def __init__(self, tolerance=1e-10, iterations=5, MLOptions={}, testUnsupported = False):
        """
        For detailed information on the possible parameters for ML, see
        http://trilinos.sandia.gov/packages/ml/documentation.html

        Currently, passing options to Aztec through ML is not supported.

        Parameters
        ----------
        tolerance : float
            Required error tolerance.
        iterations : int
            Maximum number of iterative steps to perform.
        MLOptions : dict
            Options to pass to ML.
            This will be passed to `ML.SetParameterList`.
        testUnsupported : bool
            Test smoothers that are not currently implemented in
            preconditioner objects.
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
