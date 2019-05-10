from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

__all__ = []

from fipy.terms.sourceTerm import SourceTerm

class _ExplicitSourceTerm(SourceTerm):
    r"""

    The `_ExplicitSourceTerm` discretization is given by

    .. math::

       \int_V S \,dV \simeq S_P V_P

    where :math:`S` is the `coeff` value. This source is added to the RHS vector and
    does not contribute to the solution matrix.

    """

    def _getWeight(self, var, transientGeomCoeff=None, diffusionGeomCoeff=None):
        return {
            'b vector': -1,
            'new value': 0,
            'old value': 0,
            'diagonal' : 0
        }

    def __repr__(self):
        return repr(self.coeff)
