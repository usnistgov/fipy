from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy.terms.abstractDiffusionTerm import _AbstractDiffusionTerm

__all__ = ["ExplicitDiffusionTerm"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class ExplicitDiffusionTerm(_AbstractDiffusionTerm):
    r"""
    The discretization for the `ExplicitDiffusionTerm` is given by

    .. math::

       \int_V \nabla \cdot (\Gamma\nabla\phi) dV \simeq \sum_f \Gamma_f
       \frac{\phi_A^\text{old}-\phi_P^\text{old}}{d_{AP}} A_f

    where :math:`\phi_A^\text{old}` and :math:`\phi_P^\text{old}` are the old values of the
    variable. The term is added to the RHS vector and makes no contribution to
    the solution matrix.

    """

    def _buildMatrix(self, var, SparseMatrix, boundaryConditions = (), dt=None, transientGeomCoeff=None, diffusionGeomCoeff=None):
        if hasattr(var, 'old'):
            varOld = var.old
        else:
            varOld = var

        varOld, L, b = _AbstractDiffusionTerm._buildMatrix(self, varOld, SparseMatrix, boundaryConditions = boundaryConditions, dt = dt,
                                                  transientGeomCoeff=transientGeomCoeff, diffusionGeomCoeff=diffusionGeomCoeff)

        return (var, SparseMatrix(mesh=var.mesh), b - L * var.value)

    def _getNormals(self, mesh):
        return mesh._faceCellToCellNormals

    def _treatMeshAsOrthogonal(self, mesh):
        return mesh._isOrthogonal
