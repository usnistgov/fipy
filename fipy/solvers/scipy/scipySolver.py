__docformat__ = 'restructuredtext'

__all__ = ["ScipySolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

from scipy.sparse import linalg

from fipy.matrices.scipyMatrix import _ScipyMeshMatrix
from fipy.solvers.solver import Solver
from fipy.tools import numerix

class ScipySolver(Solver):
    """
    The base `ScipySolver` class.

    .. attention:: This class is abstract. Always create one of its subclasses.
    """

    def __init__(self, tolerance="default", absolute_tolerance=0.,
                 criterion="default",
                 iterations="default", precon="default"):
        """
        Create a `Solver` object.

        Parameters
        ----------
        tolerance : float
            Required relative error tolerance.
        absolute_tolerance : float
            Required absolute error tolerance.
        criterion : {'default', 'unscaled', 'RHS', 'matrix', 'initial', 'legacy', }
            Interpretation of ``tolerance``.
            See :ref:`CONVERGENCE` for more information.
        iterations : int
            Maximum number of iterative steps to perform.
        precon
            Preconditioner to use.
        """
        self.absolute_tolerance = absolute_tolerance
        super(ScipySolver, self).__init__(tolerance=tolerance, criterion=criterion,
                                          iterations=iterations, precon=precon)

    @property
    def _matrixClass(self):
        return _ScipyMeshMatrix

    def _rhsNorm(self, L, x, b):
        return numerix.L2norm(b)

    def _matrixNorm(self, L, x, b):
        return linalg.norm(L, ord=numerix.inf)

    def _residualVectorAndNorm(self, L, x, b):
        residualVector = L * x - b

        return residualVector, numerix.L2norm(residualVector)
