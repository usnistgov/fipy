from __future__ import division
from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy.tools import numerix
from fipy.terms.asymmetricConvectionTerm import _AsymmetricConvectionTerm
from fipy.variables.faceVariable import FaceVariable

__all__ = ["HybridConvectionTerm"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class _HybridConvectionTermAlpha(FaceVariable):
    def __init__(self, P):
        FaceVariable.__init__(self, P.mesh)
        self.P = self._requires(P)

    def _calcValue(self):
        eps = 1e-3
        P  = self.P

        alpha = numerix.where(                                 P > 2., (P - 1) / P,    0.)
        alpha = numerix.where( numerix.logical_and(2. >= P, P >= -2.),         0.5, alpha)
        alpha = numerix.where(                               -2. >  P,      -1 / P, alpha)

        return alpha

class HybridConvectionTerm(_AsymmetricConvectionTerm):
    r"""
    The discretization for this :class:`~fipy.terms.term.Term` is given by

    .. math::

       \int_V \nabla \cdot (\vec{u} \phi)\,dV \simeq \sum_{f} (\vec{n}
       \cdot \vec{u})_f \phi_f A_f

    where :math:`\phi_f=\alpha_f \phi_P +(1-\alpha_f)\phi_A` and
    :math:`\alpha_f` is calculated using the hybrid scheme.
    For further details see :ref:`sec:NumericalSchemes`.
    """
    def _alpha(self, P):
        return _HybridConvectionTermAlpha(P)
