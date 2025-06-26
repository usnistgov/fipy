from __future__ import unicode_literals

from fipy.solvers.preconditioner import SolverModifyingPreconditioner

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
           "PyAMGXPreconditioner"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class PyAMGXPreconditioner(SolverModifyingPreconditioner):
    """Interface to pyamgx_ `preconditioner configuration`_ for :class:`~fipy.solvers.pyamgx.pyAMGXSolver.PyAMGXSolver`.

    .. _pyamgx: https://pyamgx.readthedocs.io
    .. _preconditioner configuration: https://pyamgx.readthedocs.io/en/latest/basic.html#config-objects
    """

    def __init__(self, **kwargs):
        """
        Parameters
        ----------
        **kwargs : dict, optional
            Extra arguments to preconditioner: refer to `preconditioner
            configuration`_ for information about possible arguments.
        """
        self.config_dict = {
            "solver": self.pctype,
            "max_iters": 1
        }
        self.config_dict.update(kwargs)

    def _applyToSolver(self, solver, matrix=None):
        solver["preconditioner"] = self.config_dict.copy()

class AMGPreconditioner(PyAMGXPreconditioner):
    """Adaptive Multigrid preconditioner for :class:`~fipy.solvers.pyamgx.pyAMGXSolver.PyAMGXSolver`.
    """

    pctype = "AMG"

class AggregationAMGPreconditioner(AMGPreconditioner):
    """Aggregation Adaptive Multigrid preconditioner for :class:`~fipy.solvers.pyamgx.pyAMGXSolver.PyAMGXSolver`.
    """

    def __init__(self):
        super(ClassicalAMGPreconditioner, self).__init__(algorithm="AGGREGATION",
                                                         selector="SIZE_2")
class BiCGStabPreconditioner(PyAMGXPreconditioner):
    """Biconjugate Gradient Stabilized preconditioner for :class:`~fipy.solvers.pyamgx.pyAMGXSolver.PyAMGXSolver`.
    """

    pctype = "PCIBCGSTAB"

class CGPreconditioner(PyAMGXPreconditioner):
    """Conjugate Gradient preconditioner for :class:`~fipy.solvers.pyamgx.pyAMGXSolver.PyAMGXSolver`.
    """

    pctype = "PCG"

class DILUPreconditioner(PyAMGXPreconditioner):
    """Diagonal Incomplete LU preconditioner for :class:`~fipy.solvers.pyamgx.pyAMGXSolver.PyAMGXSolver`.
    """

    pctype = "MULTICOLOR_DILU"

class FGMRESPreconditioner(PyAMGXPreconditioner):
    """Flexible Generalized Minimum Residual preconditioner for :class:`~fipy.solvers.pyamgx.pyAMGXSolver.PyAMGXSolver`.
    """

    pctype = "FGMRES"

class GaussSeidelPreconditioner(PyAMGXPreconditioner):
    """Gauss-Seidel preconditioner for :class:`~fipy.solvers.pyamgx.pyAMGXSolver.PyAMGXSolver`.
    """

    pctype = "MULTICOLOR_GS"

class ILUPreconditioner(PyAMGXPreconditioner):
    """Incomplete LU preconditioner for :class:`~fipy.solvers.pyamgx.pyAMGXSolver.PyAMGXSolver`.
    """

    pctype = "MULTICOLOR_GS"

class JacobiPreconditioner(PyAMGXPreconditioner):
    """Block Jacobi preconditioner for :class:`~fipy.solvers.pyamgx.pyAMGXSolver.PyAMGXSolver`.
    """

    pctype = "BLOCK_JACOBI"

class PolynomialPreconditioner(PyAMGXPreconditioner):
    """Polynomial preconditioner  for :class:`~fipy.solvers.pyamgx.pyAMGXSolver.PyAMGXSolver`.
    """

    pctype = "POLYNOMIAL"
