from __future__ import unicode_literals

from fipy.solvers.preconditioner import MatrixModifyingPreconditioner

__all__ = ["PysparsePreconditioner"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class PysparsePreconditioner(MatrixModifyingPreconditioner):
    """Base class for preconditioners of :class:`~fipy.solvers.pysparse.pysparseSolver.PysparseSolver`.

    .. attention:: This class is abstract. Always create one of its subclasses.
    """

    pass
