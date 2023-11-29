"""Boundary condition of order 0
"""
from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

"""Fixed value (Dirichlet) boundary condition
"""
__docformat__ = 'restructuredtext'

from fipy.tools import numerix

from fipy.boundaryConditions.boundaryCondition import BoundaryCondition
from fipy.tools import numerix
from fipy.tools import vector
from fipy.variables.variable import Variable

__all__ = ["FixedValue"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class FixedValue(BoundaryCondition):
    r"""Adds a Dirichlet contribution to the system of equations.

    Implements

    .. math::

       \phi|_\text{faces} = \text{value}

    The contributions are given by
    :math:`-\mathtt{value}\times G_{\text{face}}`
    for the RHS vector and :math:`G_{\text{face}}` for
    the coefficient matrix. The parameter :math:`G_{\text{face}}` represents the
    term's geometric coefficient, which depends on the type of term and the mesh
    geometry.

    Contributions are only added to entries corresponding to the
    specified faces.
    """


    def _buildMatrix(self, SparseMatrix, Ncells, MaxFaces, coeff):
        """Set boundary equal to value.

        A `tuple` of (`LL`, `bb`) is calculated, to be added to the
        Term's (:math:`\mathsf{L}`, :math:`\mathsf{b}`) matrices.

        Parameters
        ----------
        SparseMatrix : ~fipy.matrices.sparseMatrix._SparseMatrix
        Ncells : int
            Size of matrices
        MaxFaces : int
            Maximum number of faces per cell (determines number of
            non-zeros per row of :math:`\mathsf{L}`)
        coeff : list
            Contribution to adjacent cell diagonal and
            :math:`\mathsf{b}` vector by this exterior face
        """
        faces = self.faces.value

        LL = SparseMatrix(mesh=self.faces.mesh, nonZerosPerRow=1)
        LL.addAt(coeff['cell 1 diag'][faces], self.adjacentCellIDs, self.adjacentCellIDs)

        ## The following has been commented out because
        ## FixedValue's _buildMatrix() method is called for
        ## each term in the equation. Thus minusCoeff can be different for each term.
        ##
        ## if not hasattr(self, 'minusCoeff'):
        ##     self.minusCoeff = -coeff['cell 1 offdiag']
        ##     self.minusCoeff.dontCacheMe()

        bb = numerix.zeros((Ncells,), 'd')

        value = self.value
        if isinstance(value, Variable):
            value = value.value
        if value.shape == faces.shape:
            value = value[faces]

        vector.putAdd(bb, self.adjacentCellIDs, -coeff['cell 1 offdiag'].value[faces] * value)

        return (LL, bb)
