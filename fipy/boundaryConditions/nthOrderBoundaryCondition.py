__docformat__ = 'restructuredtext'

from fipy.tools import numerix

from fipy.boundaryConditions.boundaryCondition import BoundaryCondition
from fipy.boundaryConditions.fixedFlux import FixedFlux
from fipy.boundaryConditions.fixedValue import FixedValue

__all__ = ["NthOrderBoundaryCondition"]

class NthOrderBoundaryCondition(BoundaryCondition):
    """

    This boundary condition is generally used in conjunction with a
    `ImplicitDiffusionTerm` that has multiple coefficients.  It does not
    have any direct effect on the solution matrices, but its derivatives
    do.

    """

    def __init__(self, faces, value, order):
        """
        Creates an `NthOrderBoundaryCondition`.

        :Parameters:
          - `faces`: A `list` or `tuple` of `Face` objects to which this condition applies.
          - `value`: The value to impose.
          - `order`: The order of the boundary condition. An `order` of `0`
            corresponds to a `FixedValue` and an `order` of `1` corresponds to
            a `FixedFlux`. Even and odd orders behave like `FixedValue` and `FixedFlux` objects,
            respectively, but apply to higher order terms.


        """
        self.order = order
        self.derivative = {}
        BoundaryCondition.__init__(self,faces,value)

    def _buildMatrix(self, SparseMatrix, Ncells, MaxFaces, coeff):
        """Leave **L** and **b** unchanged

        :Parameters:
          - `SparseMatrix`: *unused*
          - `Ncells`:       *unused*
          - `MaxFaces`:     *unused*
          - `coeff`:        *unused*
        """
        return (0, 0)

    def _getDerivative(self, order):
        newOrder = self.order - order
        if newOrder not in self.derivative:
            if newOrder > 1:
                self.derivative[newOrder] = NthOrderBoundaryCondition(self.faces, self.value, newOrder)
            elif newOrder == 1:
                self.derivative[newOrder] = FixedFlux(self.faces, self.value)
            elif newOrder == 0:
                self.derivative[newOrder] = FixedValue(self.faces, self.value)
            else:
                self.derivative[newOrder] = None

        return self.derivative[newOrder]
