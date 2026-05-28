__docformat__ = 'restructuredtext'

from fipy.terms.explicitSourceTerm import _ExplicitSourceTerm
from fipy.terms import TermMultiplyError
from fipy.variables.cellVariable import CellVariable

__all__ = ["ResidualTerm"]

class ResidualTerm(_ExplicitSourceTerm):
    r"""

    The `ResidualTerm` is a special form of explicit `SourceTerm` that adds the
    residual of one equation to another equation. Useful for Newton's method.

    A `ResidualTerm` can be negated or multiplied by a scalar; the residual
    contribution is signed/scaled accordingly and the equation form is left
    untouched (issue #885).  This makes patterns like

       >>> import numpy as np
       >>> from fipy import Grid1D, CellVariable, TransientTerm, DiffusionTerm
       >>> from fipy.terms.residualTerm import ResidualTerm
       >>> from fipy.matrices.scipyMatrix import _ScipyMeshMatrix
       >>> mesh = Grid1D(nx=10, dx=1.)
       >>> v = CellVariable(mesh=mesh, hasOld=True,
       ...                  value=np.arange(10, dtype=float))
       >>> v.updateOld()
       >>> eq = TransientTerm(var=v) == DiffusionTerm(var=v, coeff=1.)
       >>> rt = ResidualTerm(equation=eq)

    valid for use with subtraction or scalar pre-multiplication.  Negation
    flips the sign of the residual contribution:

       >>> sandbox = CellVariable(mesh=mesh, hasOld=True, value=0.)
       >>> sandbox.updateOld()
       >>> _, _, b_pos = rt._buildAndAddMatrices(sandbox, _ScipyMeshMatrix, dt=0.1)
       >>> _, _, b_neg = (-rt)._buildAndAddMatrices(sandbox, _ScipyMeshMatrix, dt=0.1)
       >>> bool(np.allclose(np.asarray(b_pos), -np.asarray(b_neg)))
       True

    and a scalar multiplier scales it linearly:

       >>> _, _, b_3 = (3. * rt)._buildAndAddMatrices(sandbox, _ScipyMeshMatrix, dt=0.1)
       >>> bool(np.allclose(np.asarray(b_3), 3 * np.asarray(b_pos)))
       True
    """
    def __init__(self, equation, underRelaxation=1.):
        self.equation = equation
        self.underRelaxation = underRelaxation
        # Signed scale picked up by ``__neg__`` and ``__mul__``.  The default
        # value of 1 leaves the residual contribution unchanged; ``-RT`` or
        # ``alpha * RT`` flip / scale it without rebuilding the underlying
        # ``equation``.  See issue #885.
        self._scale = 1.

        _ExplicitSourceTerm.__init__(self, var=None)

    def __repr__(self):
        prefix = ""
        if self._scale == -1.:
            prefix = "-"
        elif self._scale != 1.:
            prefix = "{0!r} * ".format(self._scale)
        return prefix + r"$\Delta$[" + repr(self.equation) + "]"

    def __neg__(self):
        r"""Negate a :class:`ResidualTerm`.

        ``ResidualTerm`` carries an ``equation`` rather than the
        ``(coeff, var)`` pair that the generic
        :meth:`_NonDiffusionTerm.__neg__` assumes, so the inherited
        implementation cannot be used here without raising a
        :exc:`TypeError`.  Negation is therefore performed by flipping the
        sign of an internal ``_scale`` multiplier that
        :meth:`_buildMatrix` applies to the residual vector.

        Regression test for issue #885:

           >>> from fipy import Grid1D, CellVariable, TransientTerm
           >>> mesh = Grid1D(nx=2)
           >>> v = CellVariable(mesh=mesh, hasOld=True, value=1.)
           >>> eq = TransientTerm(var=v) == 0
           >>> from fipy.terms.residualTerm import ResidualTerm
           >>> (-ResidualTerm(equation=eq))._scale
           -1.0
           >>> (-(-ResidualTerm(equation=eq)))._scale
           1.0
        """
        ret = self.__class__(equation=self.equation,
                             underRelaxation=self.underRelaxation)
        ret._scale = -self._scale
        return ret

    def __mul__(self, other):
        r"""Scale a :class:`ResidualTerm` by a number.

        As with :meth:`__neg__`, the inherited
        :meth:`_NonDiffusionTerm.__mul__` cannot be used because
        ``ResidualTerm`` does not accept ``coeff`` in its constructor.
        Scaling is folded into the internal ``_scale`` multiplier.

        Regression test for issue #885:

           >>> from fipy import Grid1D, CellVariable, TransientTerm
           >>> mesh = Grid1D(nx=2)
           >>> v = CellVariable(mesh=mesh, hasOld=True, value=1.)
           >>> eq = TransientTerm(var=v) == 0
           >>> from fipy.terms.residualTerm import ResidualTerm
           >>> (2.5 * ResidualTerm(equation=eq))._scale
           2.5
        """
        if isinstance(other, (int, float)):
            ret = self.__class__(equation=self.equation,
                                 underRelaxation=self.underRelaxation)
            ret._scale = other * self._scale
            return ret
        else:
            raise TermMultiplyError

    __rmul__ = __mul__

    def _getGeomCoeff(self, var):
        return self.coeff

    def _buildMatrix(self, var, SparseMatrix, boundaryConditions=(), dt=None, transientGeomCoeff=None, diffusionGeomCoeff=None):
        vec = self.equation.justResidualVector(var=None,
                                               boundaryConditions=boundaryConditions,
                                               dt=dt)

        self.coeff = CellVariable(mesh=var.mesh,
                                  value=vec * self.underRelaxation * self._scale)
        self.geomCoeff = None
        self.coeffVectors = None

        return _ExplicitSourceTerm._buildMatrix(self, var=var, SparseMatrix=SparseMatrix, boundaryConditions=boundaryConditions, dt=dt, transientGeomCoeff=transientGeomCoeff, diffusionGeomCoeff=diffusionGeomCoeff)


def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()


if __name__ == "__main__":
    _test()
