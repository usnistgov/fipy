from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy.terms.abstractUpwindConvectionTerm import _AbstractUpwindConvectionTerm

__all__ = ["UpwindConvectionTerm"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class UpwindConvectionTerm(_AbstractUpwindConvectionTerm):
    r"""
    The discretization for this :class:`~fipy.terms.term.Term` is given by

    .. math::

       \int_V \nabla \cdot (\vec{u} \phi)\,dV \simeq \sum_{f} (\vec{n}
       \cdot \vec{u})_f \phi_f A_f

    where :math:`\phi_f=\alpha_f \phi_P +(1-\alpha_f)\phi_A` and
    :math:`\alpha_f` is calculated using the upwind convection scheme.
    For further details see :ref:`sec:NumericalSchemes`.
    """

    def _getDefaultSolver(self, var, solver, *args, **kwargs):
        solver = solver or super(UpwindConvectionTerm, self)._getDefaultSolver(var, solver, *args, **kwargs)
        if solver and not solver._canSolveAsymmetric():
            import warnings
            warnings.warn("%s cannot solve asymmetric matrices" % solver)
        from fipy.solvers import DefaultAsymmetricSolver
        return solver or DefaultAsymmetricSolver(*args, **kwargs)

    def _testPecletSign(self):
        r"""
            >>> from fipy import *
            >>> nx = 11
            >>> L = 1.
            >>> mesh = Grid1D(nx=nx, dx=L / nx)
            >>> var = CellVariable(mesh=mesh)

            >>> from fipy.variables.faceVariable import FaceVariable
            >>> convCoeff = FaceVariable(mesh=mesh, rank=1)
            >>> diffCoeff = FaceVariable(mesh=mesh, value=1e-20)
            >>> x = mesh.faceCenters[0]
            >>> D = 1000.
            >>> u = -1e+6
            >>> convCoeff.setValue((u,), where=x < L / 2)
            >>> diffCoeff.setValue(D, where=x > L / 2)

            >>> var.constrain(1., mesh.facesRight)
            >>> var.constrain(10, mesh.facesLeft)

            >>> dTerm = DiffusionTerm(diffCoeff)

            >>> x = mesh.cellCenters[0]
            >>> var0 = 2 * D / (-u * L + 2 * D)
            >>> analytical = 2 * (1 - var0) * x / L + 2 * var0 - 1
            >>> analytical = var0 * (x < L / 2) + analytical * (x >= L / 2)

            >>> var[:] = 0
            >>> eqn = UpwindConvectionTerm(coeff=convCoeff) == dTerm
            >>> eqn.solve(var)
            >>> print(var.allclose(analytical))
            1

            >>> var[:] = 0
            >>> eqn = TransientTerm(1e-10) == UpwindConvectionTerm(coeff=-convCoeff) +  dTerm
            >>> eqn.solve(var, dt = 1e+10)
            >>> print(var.allclose(analytical))
            1

            >>> var[:] = 0
            >>> eqn = 0 == UpwindConvectionTerm(coeff=-convCoeff) +  dTerm
            >>> eqn.solve(var)
            >>> print(var.allclose(analytical))
            1

            >>> var[:] = 0
            >>> eqn = 0 == -UpwindConvectionTerm(coeff=convCoeff) +  dTerm
            >>> eqn.solve(var)
            >>> print(var.allclose(analytical))
            1

        """
        pass

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()

