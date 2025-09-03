__docformat__ = 'restructuredtext'

from fipy.solvers.preconditioner import MatrixModifyingPreconditioner

__all__ = ["ScipyPreconditioner"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class ScipyPreconditioner(MatrixModifyingPreconditioner):
    """Base class for preconditioners for :class:`~fipy.solvers.scipy.scipySolver.ScipySolver`.

    .. attention:: This class is abstract. Always create one of its subclasses.
    """

    pass
