from __future__ import unicode_literals
from builtins import object
import copy

__all__ = ["AMGPreconditioner",
           "AggregationAMGPreconditioner",
           "BiCGStabPreconditioner",
           "CGPreconditioner",
           "DILUPreconditioner",
           "FGMRESPreconditioner",
           "GaussSeidelPreconditioner",
           "ILUPreconditioner",
           "JacobiPreconditioner",
           "PolynomialPreconditioner",
           "Preconditioner"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class Preconditioner(object):
    """Interface to pyamgx_ `preconditioner configuration`_.

    .. _pyamgx: https://pyamgx.readthedocs.io
    .. _preconditioner configuration: https://pyamgx.readthedocs.io/en/latest/basic.html#config-objects
    """

    def __init__(self, **kwargs):
        self.config_dict = {
            "solver": self.pctype,
            "max_iters": 1
        }
        self.config_dict.update(kwargs)

    def __call__(self, **kwargs):
        self.config_dict.update(kwargs)
        return copy.copy(self.config_dict)

class AMGPreconditioner(Preconditioner):
    """
    Adaptive Multigrid Preconditioner for pyamgx solvers.

    """
    pctype = "AMG"

class AggregationAMGPreconditioner(AMGPreconditioner):
    """
    Aggregation Adaptive Multigrid Preconditioner for pyamgx solvers.

    """
    def __init__(self):
        super(ClassicalAMGPreconditioner, self).__init__(algorithm="AGGREGATION",
                                                         selector="SIZE_2")
class BiCGStabPreconditioner(Preconditioner):
    """
    Biconjugate Gradient Stabilized Preconditioner for pyamgx solvers.

    """
    pctype = "PCIBCGSTAB"

class CGPreconditioner(Preconditioner):
    """
    Conjugate Gradient Preconditioner for pyamgx solvers.

    """
    pctype = "PCG"

class DILUPreconditioner(Preconditioner):
    """
    DILU Preconditioner for pyamgx solvers.

    """
    pctype = "MULTICOLOR_DILU"

class FGMRESPreconditioner(Preconditioner):
    """
    Flexible Generalized Mimumal Residual Preconditioner for pyamgx solvers.

    """
    pctype = "FGMRES"

class GaussSeidelPreconditioner(Preconditioner):
    """
    Gauss-Seidel Preconditioner for pyamgx solvers.

    """
    pctype = "MULTICOLOR_GS"

class ILUPreconditioner(Preconditioner):
    """
    ILU Preconditioner for pyamgx solvers.

    """
    pctype = "MULTICOLOR_GS"

class JacobiPreconditioner(Preconditioner):
    """
    Block Jacobi Preconditioner for pyamgx solvers.

    """
    pctype = "BLOCK_JACOBI"

class PolynomialPreconditioner(Preconditioner):
    """
    Polynomial Preconditioner for pyamgx solvers.

    """
    pctype = "POLYNOMIAL"
