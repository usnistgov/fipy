from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy.terms.sourceTerm import SourceTerm
from fipy.tools import numerix

__all__ = ["ImplicitSourceTerm"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class ImplicitSourceTerm(SourceTerm):
    r"""

    The `ImplicitSourceTerm` represents

    .. math::

       \int_V \phi S \,dV \simeq \phi_P S_P V_P

    where :math:`S` is the `coeff` value.
    """

    def __init__(self, coeff=1., var=None):
        r"""
        Parameters
        ----------
        coeff : float or ~fipy.variables.cellVariable.CellVariable
            Proportionality coefficient :math:`S` (default: 1)
        var : ~fipy.variables.cellVariable.CellVariable
            Variable :math:`\phi` that
            :class:`~fipy.terms.implicitSourceTerm.ImplicitSourceTerm` is
            implicit in.
        """
        super(ImplicitSourceTerm, self).__init__(coeff=coeff, var=var)

    def _getWeight(self, var, transientGeomCoeff=None, diffusionGeomCoeff=None):
        """
        Test for a bug due to the sign operator not being updating
        correctly.

            >>> from fipy import *
            >>> m = Grid1D(nx=1)
            >>> v = CellVariable(mesh=m, value=1.)
            >>> eq = TransientTerm() == ImplicitSourceTerm(v)
            >>> eq.solve(v, dt=1.)
            >>> print(v)
            [ 2.]
            >>> v.setValue(-1.)
            >>> eq.solve(v, dt=1.)
            >>> print(v)
            [-0.5]

        """

        coeff = self._getGeomCoeff(var)
        diagonalSign = self._getDiagonalSign(transientGeomCoeff, diffusionGeomCoeff)
        combinedSign = numerix.array(diagonalSign)[..., numerix.newaxis] * numerix.sign(coeff)

        return {'diagonal' : (combinedSign >= 0),
                'old value' : numerix.zeros(var.shape, 'd'),
                'b vector' :  -var * (combinedSign < 0),
                'new value' : numerix.zeros(var.shape, 'd')}

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()

