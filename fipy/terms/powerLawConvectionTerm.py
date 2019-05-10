from __future__ import division
from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy.terms.asymmetricConvectionTerm import _AsymmetricConvectionTerm
from fipy.variables.faceVariable import FaceVariable
from fipy.tools.dimensions.physicalField import PhysicalField
from fipy.tools import inline
from fipy.tools import numerix

__all__ = ["PowerLawConvectionTerm"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class _PowerLawConvectionTermAlpha(FaceVariable):
    """

    Test case added because `and` was being used instead of bitwise `&`.

        >>> from fipy.meshes import Grid1D
        >>> mesh = Grid1D(nx = 3)
        >>> from fipy.variables.faceVariable import FaceVariable
        >>> P = FaceVariable(mesh=mesh, value=(1e-3, 1e+71, 1e-3, 1e-3))
        >>>
        >>> alpha = PowerLawConvectionTerm([1])._alpha(P)
        >>> print(numerix.allclose(alpha, [ 0.5,  1.,   0.5, 0.5]))
        True
    """
    def __init__(self, P):
        FaceVariable.__init__(self, mesh = P.mesh, elementshape=P.shape[:-1])
        self.P = self._requires(P)
        self.eps = 1e-3

    if inline.doInline:
        def _calcValue(self):
            eps = self.eps
            P = self.P.numericValue
            alpha = self._array.copy()

            inline._runInline("""
                if (fabs(P[i]) < eps) {
                    P[i] = eps;
                }

                alpha[i] = 0.5;

                if (P[i] > 10.) {
                    alpha[i] = (P[i] - 1.) / P[i];
                } else if (10. >= P[i] && P[i] > eps) {
                    double	tmp = (1. - P[i] / 10.);
                    double	tmpSqr = tmp * tmp;
                    alpha[i] = ((P[i] - 1.) + tmpSqr*tmpSqr*tmp) / P[i];
                } else if (-eps > P[i] && P[i] >= -10.) {
                    double	tmp = (1. + P[i] / 10.);
                    double	tmpSqr = tmp * tmp;
                    alpha[i] = (tmpSqr*tmpSqr*tmp - 1.) / P[i];
                } else if (P[i] < -10.) {
                    alpha[i] = -1. / P[i];
                }
            """,
            alpha = alpha, eps = eps, P = P,
            ni = len(P.flat)
            )

            return self._makeValue(value = alpha)
    else:
        def _calcValue(self):
            eps = self.eps
            P = self.P.numericValue
            P = numerix.where(abs(P) < eps, eps, P)

            alpha = numerix.where(                  P > 10.,                     (P - 1.) / P,   0.5)

            tmp = (1. - P / 10.)
            tmpSqr = tmp * tmp
            alpha = numerix.where(  (10. >= P) & (P > eps),  ((P-1.) + tmpSqr*tmpSqr*tmp) / P, alpha)

            tmp = (1. + P / 10.)
            tmpSqr = tmp * tmp
            alpha = numerix.where((-eps >  P) & (P >= -10.),     (tmpSqr*tmpSqr*tmp - 1.) / P, alpha)

            alpha = numerix.where(                 P < -10.,                          -1. / P, alpha)

            return PhysicalField(value = alpha)

class PowerLawConvectionTerm(_AsymmetricConvectionTerm):
    r"""
    The discretization for this :class:`~fipy.terms.term.Term` is given by

    .. math::

       \int_V \nabla \cdot (\vec{u} \phi)\,dV \simeq \sum_{f} (\vec{n}
       \cdot \vec{u})_f \phi_f A_f

    where :math:`\phi_f=\alpha_f \phi_P +(1-\alpha_f)\phi_A` and
    :math:`\alpha_f` is calculated using the power law scheme.
    For further details see :ref:`sec:NumericalSchemes`.
    """
    def _alpha(self, P):
        return _PowerLawConvectionTermAlpha(P)

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()


