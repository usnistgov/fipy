__docformat__ = 'restructuredtext'

from fipy.tools import numerix

from fipy.variables.cellToFaceVariable import _CellToFaceVariable
from fipy.tools.numerix import exp, expm1, minimum, take, where

__all__ = ["ScharfetterGummelFaceVariable"]

class ScharfetterGummelFaceVariable(_CellToFaceVariable):
    r"""Scharfetter-Gummel interpolation of :math:`e^\phi` to the faces.

    For the cell values :math:`a` and :math:`b` on either side of a face,
    the face value is

    .. math::

       \frac{(a - b) e^a}{e^{a - b} - 1}
       = e^{\min(a, b)} \frac{d}{1 - e^{-d}}, \qquad d = |a - b|,

    which tends to :math:`e^a` as :math:`a - b \to 0`. The right-hand form
    is symmetric in :math:`a` and :math:`b` and is evaluated with
    :func:`expm1`, so a large difference in either direction gives the
    finite value.

    >>> from fipy import CellVariable, Grid1D
    >>> mesh = Grid1D(nx=2, dx=1.)
    >>> for values in [(100., 0.), (0., 100.)]:
    ...     phi = CellVariable(mesh=mesh, value=values)
    ...     print(numerix.allclose(ScharfetterGummelFaceVariable(phi)[1],
    ...                            100. / -numerix.expm1(-100.)))
    True
    True

    Equal values give :math:`e^a`, and differences just below, at and above
    100 agree in both orientations with the closed form.

    >>> phi = CellVariable(mesh=mesh, value=(3., 3.))
    >>> print(numerix.allclose(ScharfetterGummelFaceVariable(phi)[1],
    ...                        numerix.exp(3.)))
    True
    >>> for d in [1e-20, 1e-8, 0.5, 99.9, 100., 100.1, 700.]:
    ...     expected = d / -numerix.expm1(-d)
    ...     for values in [(d, 0.), (0., d)]:
    ...         phi = CellVariable(mesh=mesh, value=values)
    ...         value = ScharfetterGummelFaceVariable(phi)[1]
    ...         assert numerix.allclose(value, expected, rtol=1e-12), (d, values)

    Reflecting a diffusion problem with this coefficient reflects its flux.

    >>> from fipy import DiffusionTerm
    >>> fluxes = []
    >>> for values, left, right in [((100., 0.), 0., 1.), ((0., 100.), 1., 0.)]:
    ...     phi = CellVariable(mesh=mesh, value=values)
    ...     coeff = ScharfetterGummelFaceVariable(phi)
    ...     u = CellVariable(mesh=mesh, value=0.)
    ...     u.constrain(left, mesh.facesLeft)
    ...     u.constrain(right, mesh.facesRight)
    ...     DiffusionTerm(coeff=coeff).solve(var=u)
    ...     fluxes.append(abs(coeff.value[1] * (u.value[1] - u.value[0])))
    >>> print(numerix.allclose(fluxes[0], fluxes[1]))
    True
    """

    def __init__(self, var, boundaryConditions=()):
        _CellToFaceVariable.__init__(self, var)
        self.bcs = boundaryConditions

    def _calcValue_(self, alpha, id1, id2):
        cell1 = take(self.var, id1)
        cell2 = take(self.var, id2)
        delta = abs(cell1 - cell2)

        safe = where(delta == 0, 1., delta)
        value = where(delta == 0, 1., safe / -expm1(-safe))
        value = value * exp(minimum(cell1, cell2))

        for bc in self.bcs:
            if isinstance(bc, FixedValue):
                value[bc.faces.value] = bc._value

        return value

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
