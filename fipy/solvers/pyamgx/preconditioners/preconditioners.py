from __future__ import unicode_literals
from builtins import object
import copy

__all__ = ["AggregationAMGPreconditioner",
           "ClassicalAMGPreconditioner",
           "CGPreconditioner",
           "BiCGStabPreconditioner",
           "FGMRESPreconditioner",
           "BlockJacobiPreconditioner",
           "MultiColorDILUPreconditioner",
           "PolynomialPreconditioner",
           "MultiColorGSPreconditioner",
           "Preconditioner"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class Preconditioner(object):
    """Interface to pyamgx_ `preconditioner configuration`_.

    .. _pyamgx: https://pyamgx.readthedocs.io
    .. _preconditioner configuration: https://pyamgx.readthedocs.io/en/latest/basic.html#config-objects
    """

    def __init__(self, preconditioner_type, **kwargs):
        self.config_dict = {
            "solver": preconditioner_type,
            "max_iters": 1
        }
        self.config_dict.update(kwargs)
    def __call__(self, **kwargs):
        """
        Parameters
        ----------
        **kwargs
            Other AMGX solver options
        """
        self.config_dict.update(kwargs)
        return copy.copy(self.config_dict)

AggregationAMGPreconditioner = Preconditioner("AMG", algorithm="AGGREGATION", selector="SIZE_2")
ClassicalAMGPreconditioner = Preconditioner("AMG")
CGPreconditioner = Preconditioner("PCG")
BiCGStabPreconditioner = Preconditioner("PCIBCGSTAB")
FGMRESPreconditioner = Preconditioner("FGMRES")
BlockJacobiPreconditioner = Preconditioner("BLOCK_JACOBI")
MultiColorDILUPreconditioner = Preconditioner("MULTICOLOR_DILU")
PolynomialPreconditioner = Preconditioner("POLYNOMIAL")
MultiColorGSPreconditioner = Preconditioner("MULTICOLOR_GS")
