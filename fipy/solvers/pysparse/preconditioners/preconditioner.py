from __future__ import unicode_literals
from builtins import object
__all__ = ["Preconditioner"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class Preconditioner(object):
    """
    Base preconditioner class

    .. attention:: This class is abstract. Always create one of its subclasses.
    """

    def _applyToMatrix(self, A):
        """
        Returns the function used for Pysparse
        preconditioning.
        """
        raise NotImplementedError
