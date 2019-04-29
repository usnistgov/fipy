from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy.terms.abstractUpwindConvectionTerm import _AbstractUpwindConvectionTerm
from fipy.tools import numerix
from fipy.terms import TransientTermError

__all__ = ["ExplicitUpwindConvectionTerm"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class ExplicitUpwindConvectionTerm(_AbstractUpwindConvectionTerm):
    r"""
    The discretization for this :class:`~fipy.terms.term.Term` is given by

    .. math::

       \int_V \nabla \cdot (\vec{u} \phi)\,dV \simeq \sum_{f} (\vec{n}
       \cdot \vec{u})_f \phi_f A_f

    where :math:`\phi_f=\alpha_f \phi_P^\text{old} +(1-\alpha_f)\phi_A^\text{old}` and
    :math:`\alpha_f` is calculated using the upwind scheme.
    For further details see :ref:`sec:NumericalSchemes`.
    """

    def _getOldAdjacentValues(self, oldArray, id1, id2, dt):
        if dt is None:
            raise TransientTermError
        return numerix.take(oldArray, id1), numerix.take(oldArray, id2)

    def _getWeight(self, var, transientGeomCoeff=None, diffusionGeomCoeff=None):
        weight = _AbstractUpwindConvectionTerm._getWeight(self, var, transientGeomCoeff, diffusionGeomCoeff)
        if 'implicit' in list(weight.keys()):
            weight['explicit'] = weight['implicit']
            del weight['implicit']

        return weight
