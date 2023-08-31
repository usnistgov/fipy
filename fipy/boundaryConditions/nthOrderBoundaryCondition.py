"""Boundary condition of specified derivative order
"""
from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy.tools import numerix

from fipy.boundaryConditions.boundaryCondition import BoundaryCondition
from fipy.boundaryConditions.fixedFlux import FixedFlux
from fipy.boundaryConditions.fixedValue import FixedValue

__all__ = ["NthOrderBoundaryCondition"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class NthOrderBoundaryCondition(BoundaryCondition):
    r"""Adds an appropriate contribution to the system of equations

    Implements

    .. math::

       \hat{n}\cdot\nabla^\text{order} \phi |_\text{faces} = \text{value}

    This boundary condition is generally used in conjunction with a
    `ImplicitDiffusionTerm` that has multiple coefficients.  It does not
    have any direct effect on the solution matrices, but its derivatives
    do.
    """

    def __init__(self, faces, value, order):
        """
        Creates an `NthOrderBoundaryCondition`.

        Parameters
        ----------
        faces : :obj:`~fipy.variables.faceVariable.FaceVariable` of :obj:`bool`
            Mask of faces where this condition applies.
        value : float
            Value to impose.
        order : int
            Order of the boundary condition. An `order` of `0`
            corresponds to a `FixedValue` and an `order` of `1` corresponds to
            a `FixedFlux`. Even and odd orders behave like `FixedValue` and `FixedFlux` objects,
            respectively, but apply to higher order terms.
        """
        self.order = order
        self.derivative = {}
        BoundaryCondition.__init__(self, faces, value)

    def _buildMatrix(self, SparseMatrix, Ncells, MaxFaces, coeff):
        """Leave **L** and **b** unchanged

        Parameters
        ----------
        SparseMatrix
            *unused*
        Ncells
            *unused*
        MaxFaces
            *unused*
        coeff
            *unused*
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
