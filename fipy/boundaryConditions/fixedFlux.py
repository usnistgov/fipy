"""Boundary condition of order 1
"""
from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy.tools import numerix

from fipy.boundaryConditions.boundaryCondition import BoundaryCondition
from fipy.boundaryConditions.fixedValue import FixedValue
from fipy.tools import vector

__all__ = ["FixedFlux"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class FixedFlux(BoundaryCondition):
    r"""Adds a Neumann contribution to the system of equations.

    Implements

    .. math::

       \hat{n}\cdot\vec{J}|_\text{faces} = \text{value}

    The
    contribution, given by `value`, is only added to entries corresponding to
    the specified `faces`, and is weighted by the face areas.

    """

    def __init__(self, faces, value):
        """
        Parameters
        ----------
        faces : :obj:`~fipy.variables.faceVariable.FaceVariable` of :obj:`bool`
            Mask of faces where this condition applies.
        value : float
            Value to impose.
        """
        BoundaryCondition.__init__(self, faces, value)
        ## The extra index [self.faces.value] makes self.contribution the same length as self.adjacentCellIDs
        self.contribution = (self.value * self.faces.mesh._faceAreas)[self.faces.value]

    def _buildMatrix(self, SparseMatrix, Ncells, MaxFaces, coeff):
        """Leave **L** unchanged and add gradient to **b**

        Parameters
        ----------
        SparseMatrix : ~fipy.matrices.sparseMatrix._SparseMatrix
            *unused* (empty matrix)
        Ncells : int
            Size of **b** vector
        MaxFaces
            *unused*
        coeff : list
            *unused*
        """

        bb = numerix.zeros((Ncells,), 'd')

        if not self.boundaryConditionApplied:
            vector.putAdd(bb, self.adjacentCellIDs, -self.contribution)
            self.boundaryConditionApplied = True

        return (0, bb)

    def _getDerivative(self, order):
        if order == 1:
            return FixedValue(self.faces, self.value)
        else:
            return BoundaryCondition._getDerivative(self, order)

    def _test(self):
        """
        The following tests check that `self.contributions` is the same length as
        `self.adjacentCellIDs`.

           >>> from fipy import *
           >>> m = Grid1D(nx = 10)
           >>> v = FaceVariable(mesh=m)
           >>> bc = FixedFlux(value=v.globalValue[-1], faces=m.facesRight)
           >>> len(bc.contribution) == len(bc.adjacentCellIDs)
           True
           >>> bc = FixedFlux(value=v, faces=m.facesRight)
           >>> len(bc.contribution) == len(bc.adjacentCellIDs)
           True
        """

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
