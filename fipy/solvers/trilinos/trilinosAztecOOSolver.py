from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from PyTrilinos import AztecOO

_reason = {AztecOO.AZ_normal : 'AztecOO.AZ_normal',
           AztecOO.AZ_param : 'AztecOO.AZ_param',
           AztecOO.AZ_breakdown : 'AztecOO.AZ_breakdown',
           AztecOO.AZ_loss : 'AztecOO.AZ_loss',
           AztecOO.AZ_ill_cond : 'AztecOO.AZ_ill_cond',
           AztecOO.AZ_maxits : 'AztecOO.AZ_maxits'}

from fipy.solvers.trilinos.trilinosSolver import TrilinosSolver
from fipy.solvers.trilinos.preconditioners.jacobiPreconditioner import JacobiPreconditioner
from fipy.tools.timer import Timer

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

        self._log.debug("BEGIN precondition")

        with Timer() as t:
            if self.preconditioner is not None:
                self.preconditioner._applyToSolver(solver=Solver, matrix=L)
            else:
                Solver.SetAztecOption(AztecOO.AZ_precond, AztecOO.AZ_none)

        self._log.debug("END precondition - {} ns".format(t.elapsed))

        self._log.debug("BEGIN solve")

        with Timer() as t:
            output = Solver.Iterate(self.iterations, self.tolerance)

        self._log.debug("END solve - {} ns".format(t.elapsed))

        if self.preconditioner is not None:
            if hasattr(self.preconditioner, 'Prec'):
                del self.preconditioner.Prec

        status = Solver.GetAztecStatus()
        self._log.debug('iterations: %d / %d', status[AztecOO.AZ_its], self.iterations)
        self._log.debug('failure: %s', _reason[status[AztecOO.AZ_why]])
        self._log.debug('AztecOO.AZ_r: %s', status[AztecOO.AZ_r])
        self._log.debug('AztecOO.AZ_scaled_r: %s', status[AztecOO.AZ_scaled_r])
        self._log.debug('AztecOO.AZ_solve_time: %s', status[AztecOO.AZ_solve_time])
        self._log.debug('AztecOO.AZ_Aztec_version: %s', status[AztecOO.AZ_Aztec_version])

        return output
