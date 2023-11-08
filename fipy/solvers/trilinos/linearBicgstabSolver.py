from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from PyTrilinos import AztecOO

from fipy.solvers.trilinos.trilinosAztecOOSolver import TrilinosAztecOOSolver
from fipy.solvers.trilinos.preconditioners.jacobiPreconditioner import JacobiPreconditioner

__all__ = ["LinearBicgstabSolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class LinearBicgstabSolver(TrilinosAztecOOSolver):

    """
    The `LinearBicgstabSolver` is an interface to the biconjugate gradient
    stabilized solver in Trilinos, using the `JacobiPreconditioner` by default.

    """

    def __init__(self, tolerance=1e-10, criterion="default",
                 iterations=1000, precon=JacobiPreconditioner()):
        """
        Parameters
        ----------
        tolerance : float
            Required error tolerance.
        criterion : {'default', 'initial', 'unscaled', 'RHS', 'matrix', 'solution', 'legacy'}
            Interpretation of ``tolerance``.
            See :ref:`CONVERGENCE` for more information.
        iterations : int
            Maximum number of iterative steps to perform.
        precon : ~fipy.solvers.trilinos.preconditioners.preconditioner.Preconditioner
        """
        TrilinosAztecOOSolver.__init__(self, tolerance=tolerance, criterion=criterion,
                                       iterations=iterations, precon=precon)
        self.solver = AztecOO.AZ_bicgstab
