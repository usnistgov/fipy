from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from PyTrilinos import AztecOO

from fipy.solvers.trilinos.trilinosAztecOOSolver import TrilinosAztecOOSolver
from fipy.solvers.trilinos.preconditioners.multilevelDDPreconditioner import MultilevelDDPreconditioner

__all__ = ["LinearCGSSolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class LinearCGSSolver(TrilinosAztecOOSolver):

    """
    The `LinearCGSSolver` is an interface to the CGS solver in Trilinos, using
    the `MultilevelSGSPreconditioner` by default.

    """

    solver = AztecOO.AZ_cgs

    DEFAULT_PRECONDITIONER = MultilevelDDPreconditioner
