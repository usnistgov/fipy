


__docformat__ = 'restructuredtext'

from fipy.terms.abstractConvectionTerm import _AbstractConvectionTerm
from fipy.variables.faceVariable import FaceVariable

__all__ = ["CentralDifferenceConvectionTerm"]

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
