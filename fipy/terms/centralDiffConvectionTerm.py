from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy.terms.abstractConvectionTerm import _AbstractConvectionTerm
from fipy.variables.faceVariable import FaceVariable

__all__ = ["CentralDifferenceConvectionTerm"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

from fipy.tools import numerix

class _CentralDifferenceConvectionTermAlpha(FaceVariable):
    def __init__(self, P):
        FaceVariable.__init__(self, P.mesh, value=0.5)

class CentralDifferenceConvectionTerm(_AbstractConvectionTerm):
    r"""

    This :class:`~fipy.terms.term.Term` represents

    .. math::

       \int_V \nabla \cdot (\vec{u} \phi)\,dV \simeq \sum_{f} (\vec{n}
       \cdot \vec{u})_f \phi_f A_f

    where :math:`\phi_f=\alpha_f \phi_P +(1-\alpha_f)\phi_A` and
    :math:`\alpha_f` is calculated using the central differencing scheme.
    For further details see :ref:`sec:NumericalSchemes`.
    """
    def _alpha(self, P):
        return _CentralDifferenceConvectionTermAlpha(P)

    def _buildMatrix(self, var, SparseMatrix, boundaryConditions=(), dt=None, transientGeomCoeff=None, diffusionGeomCoeff=None):

        var, L, b = _AbstractConvectionTerm._buildMatrix(self, var, SparseMatrix, boundaryConditions=boundaryConditions, dt=dt, transientGeomCoeff=transientGeomCoeff, diffusionGeomCoeff=diffusionGeomCoeff)

##        if var.rank != 1:

        mesh = var.mesh
        
        ids = self._reshapeIDs(var, numerix.arange(mesh.numberOfCells))
        if (hasattr(self, 'constraintL')) or (hasattr(self, 'constraintB')):
            # remove L from matrix and b from right hand side
            # necessary due to call of _buildMatrix on super class
            L.addAt(numerix.array(-self.constraintL).ravel(), ids.ravel(), ids.swapaxes(0, 1).ravel())
            b -= numerix.reshape(self.constraintB.value, ids.shape).sum(0).ravel() 
            
            constraintMaskFixedGradient = var.faceGrad.constraintMask
            constraintMaskFixedValue = var.arithmeticFaceValue.constraintMask
            
            weight = self._getWeight(var, transientGeomCoeff, diffusionGeomCoeff)

            if 'implicit' in weight:
                alpha = weight['implicit']['cell 1 diag']
            else:
                alpha = 0.0

            exteriorCoeff =  self.coeff * mesh.exteriorFaces

            self.constraintL = (constraintMaskFixedGradient * exteriorCoeff).divergence * mesh.cellVolumes
            self.constraintB =  -(var.arithmeticFaceValue * constraintMaskFixedValue * exteriorCoeff).divergence * mesh.cellVolumes

        L.addAt(numerix.array(self.constraintL).ravel(), ids.ravel(), ids.swapaxes(0, 1).ravel())
        b += numerix.reshape(self.constraintB.value, ids.shape).sum(0).ravel()

        return (var, L, b)
