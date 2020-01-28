from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

import os

from PyTrilinos import AztecOO

from fipy.solvers.trilinos.trilinosSolver import TrilinosSolver
from fipy.solvers.trilinos.preconditioners.jacobiPreconditioner import JacobiPreconditioner

__all__ = ["TrilinosAztecOOSolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class TrilinosAztecOOSolver(TrilinosSolver):

    """
    .. attention:: This class is abstract, always create on of its subclasses.
       It provides the code to call all solvers from the Trilinos AztecOO package.

    """

    def __init__(self, tolerance=1e-10, iterations=1000, precon=JacobiPreconditioner()):
        """
        Parameters
        ----------
        tolerance : float
            Required error tolerance.
        iterations : int
            Maximum number of iterative steps to perform.
        precon : ~fipy.solvers.trilinos.preconditioners.preconditioner.Preconditioner
        """
        if self.__class__ is TrilinosAztecOOSolver:
            raise NotImplementedError("can't instantiate abstract base class")

        TrilinosSolver.__init__(self, tolerance=tolerance,
                                iterations=iterations, precon=None)
        self.preconditioner = precon

    def _solve_(self, L, x, b):

        Solver = AztecOO.AztecOO(L, x, b)
        Solver.SetAztecOption(AztecOO.AZ_solver, self.solver)

##        Solver.SetAztecOption(AztecOO.AZ_kspace, 30)

        Solver.SetAztecOption(AztecOO.AZ_output, AztecOO.AZ_none)

        if self.preconditioner is not None:
            self.preconditioner._applyToSolver(solver=Solver, matrix=L)
        else:
            Solver.SetAztecOption(AztecOO.AZ_precond, AztecOO.AZ_none)

        output = Solver.Iterate(self.iterations, self.tolerance)

        if self.preconditioner is not None:
            if hasattr(self.preconditioner, 'Prec'):
                del self.preconditioner.Prec

        if 'FIPY_VERBOSE_SOLVER' in os.environ:
            status = Solver.GetAztecStatus()

            from fipy.tools.debug import PRINT
            PRINT('iterations: %d / %d' % (status[AztecOO.AZ_its], self.iterations))
            failure = {AztecOO.AZ_normal : 'AztecOO.AZ_normal',
                       AztecOO.AZ_param : 'AztecOO.AZ_param',
                       AztecOO.AZ_breakdown : 'AztecOO.AZ_breakdown',
                       AztecOO.AZ_loss : 'AztecOO.AZ_loss',
                       AztecOO.AZ_ill_cond : 'AztecOO.AZ_ill_cond',
                       AztecOO.AZ_maxits : 'AztecOO.AZ_maxits'}

            PRINT('failure', failure[status[AztecOO.AZ_why]])

            PRINT('AztecOO.AZ_r:', status[AztecOO.AZ_r])
            PRINT('AztecOO.AZ_scaled_r:', status[AztecOO.AZ_scaled_r])
            PRINT('AztecOO.AZ_rec_r:', status[AztecOO.AZ_rec_r])
            PRINT('AztecOO.AZ_solve_time:', status[AztecOO.AZ_solve_time])
            PRINT('AztecOO.AZ_Aztec_version:', status[AztecOO.AZ_Aztec_version])

        return output
