__docformat__ = 'restructuredtext'

from PyTrilinos import AztecOO

from fipy.solvers.trilinos.trilinosAztecOOSolver import TrilinosAztecOOSolver
from fipy.solvers.trilinos.preconditioners.multilevelDDPreconditioner import MultilevelDDPreconditioner

__all__ = ["LinearCGSSolver"]

class LinearCGSSolver(TrilinosAztecOOSolver):

    """Interface to the conjugate gradient squared (:term:`CGS`) solver in
    :ref:`TRILINOS`

    Uses the
    :class:`~fipy.solvers.trilinos.preconditioners.multilevelDDPreconditioner.MultilevelDDPreconditioner`
    by default.
    """

    solver = AztecOO.AZ_cgs

    DEFAULT_PRECONDITIONER = MultilevelDDPreconditioner
