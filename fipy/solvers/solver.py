"""
The iterative solvers may output warnings if the solution is considered
unsatisfactory. If you are not interested in these warnings, you can invoke
python with a warning filter such as::

    $ python -Wignore::fipy.SolverConvergenceWarning myscript.py

If you are extremely concerned about your preconditioner for some reason, you
can abort whenever it has problems with::

    $ python -Werror::fipy.PreconditionerWarning myscript.py

"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from builtins import object
from builtins import str
__docformat__ = 'restructuredtext'

import logging
import warnings

from fipy.tools import numerix
from .convergence import ConvergenceBase

__all__ = ["SolverConvergenceWarning", "NormalConvergence", "MaximumIterationWarning",
           "PreconditionerWarning", "IllConditionedPreconditionerWarning",
           "PreconditionerNotPositiveDefiniteWarning", "MatrixIllConditionedWarning",
           "StagnatedSolverWarning", "ScalarQuantityOutOfRangeWarning",
           "IllegalInputOrBreakdownWarning",
           "ParameterWarning", "BreakdownWarning", "LossOfPrecisionWarning",
           "Solver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class SolverConvergenceWarning(Warning):
    def __init__(self, solver, iter, relres):
        self.solver = solver
        self.iter = iter
        self.relres = relres

    def __str__(self):
        return "%s failed. Iterations: %g. Relative error: %g" % (str(self.solver), self.iter, self.relres)

class NormalConvergence(SolverConvergenceWarning):
    def __str__(self):
        return "User requested convergence criteria is satisfied. Iterations: {0}. Relative error: {1}".format(self.iter, self.relres)

class MaximumIterationWarning(SolverConvergenceWarning):
    def __str__(self):
        return "Iterations: %g. Relative error: %g" % (self.iter, self.relres)

class PreconditionerWarning(SolverConvergenceWarning):
    pass

class IllConditionedPreconditionerWarning(PreconditionerWarning):
    def __str__(self):
        return "The system involving the preconditioner was ill-conditioned. Relative error: %g" % (self.relres)

class PreconditionerNotPositiveDefiniteWarning(PreconditionerWarning):
    def __str__(self):
        return "The preconditioning matrix does not appear to be positive definite. Relative error: %g" % (self.relres)

class MatrixIllConditionedWarning(SolverConvergenceWarning):
    def __str__(self):
        return "The matrix appears to be very ill-conditioned. Relative error: %g" % (self.relres)

class StagnatedSolverWarning(SolverConvergenceWarning):
    def __str__(self):
        return "The solver stagnated. Iterations: %g. Relative error: %g" % (self.iter, self.relres)

class ScalarQuantityOutOfRangeWarning(SolverConvergenceWarning):
    def __str__(self):
        return "A scalar quantity became too small or too large to continue computing. Iterations: %g. Relative error: %g" % (self.iter, self.relres)

class IllegalInputOrBreakdownWarning(SolverConvergenceWarning):
    def __str__(self):
        return "{0} received illegal input or had a breakdown." \
          "Iterations: {1}. Relative error: {2}".format(self.solver, self.iter, self.relres)

class ParameterWarning(SolverConvergenceWarning):
    def __str__(self):
        return "User requested option is not available for {0}.".format(self.solver)

class BreakdownWarning(SolverConvergenceWarning):
    def __str__(self):
        return "Numerical breakdown occurred. Iterations: {0}. Relative error: {1}".format(self.iter, self.relres)

class LossOfPrecisionWarning(SolverConvergenceWarning):
    def __str__(self):
        return "Numerical loss of precision occurred. Iterations: {0}. Relative error: {1}".format(self.iter, self.relres)



class Solver(object):
    """
    The base `LinearXSolver` class.

    .. attention:: This class is abstract. Always create one of its subclasses.
    """

    def __init__(self, tolerance=1e-10, criterion="default",
                 iterations=1000, precon=None):
        """
        Create a `Solver` object.

        Parameters
        ----------
        tolerance : float
            Required error tolerance.
        criterion : {'default', 'initial', 'unscaled', 'RHS', 'matrix', 'solution', 'preconditioned', 'natural'}, optional
            Interpretation of ``tolerance``.
            See :ref:`CONVERGENCE` for more information.
        iterations : int
            Maximum number of iterative steps to perform.
        precon
            Preconditioner to use.  Not all solver suites support
            preconditioners.
        """
        if self.__class__ is Solver:
            raise NotImplementedError("can't instantiate abstract base class")

        self.tolerance = tolerance
        self.criterion = criterion
        self.iterations = iterations

        self.preconditioner = precon

        self._log = logging.getLogger(self.__class__.__module__
                                      + "." + self.__class__.__name__)

    def _storeMatrix(self, var, matrix, RHSvector):
        self.var = var
        self.matrix = matrix
        self.RHSvector = RHSvector

    def _solve(self):
        raise NotImplementedError

    def _solve_(self, L, x, b):
        raise NotImplementedError

    def _setConvergence(self, suite, code, iterations, residual, actual_code=None, **kwargs):
        cls = ConvergenceBase.code_registry[(suite, code)]
        self.convergence = cls(solver=self,
                               iterations=iterations,
                               residual=residual,
                               criterion=self.criterion,
                               actual_code=actual_code,
                               **kwargs)

    def _defaultNorm(self, L, x, b):
        raise NotImplementedError

    def _unscaledNorm(self, L, x, b):
        return 1.

    def _rhsNorm(self, L, x, b):
        raise NotImplementedError

    def _matrixNorm(self, L, x, b):
        raise NotImplementedError

    def _residualVectorAndNorm(self, L, x, b):
        raise NotImplementedError

    def _residualNorm(self, L, x, b):
        _, residual = self._residualVectorAndNorm(L, x, b)

        return residual

    def _solutionNorm(self, L, x, b):
        raise NotImplementedError

    def _preconditionedNorm(self, L, x, b):
        raise NotImplementedError

    def _naturalNorm(self, L, x, b):
        raise NotImplementedError

    @property
    def _Lxb(self):
        raise NotImplementedError

    @property
    def _norms(self):
        L, x, b = self._Lxb
        return (self._matrixNorm(L, x, b),
                self._rhsNorm(L, x, b),
                self._residualNorm(L, x, b))

    def _adaptDefaultTolerance(self, L, x, b):
        raise NotImplementedError

    def _adaptUnscaledTolerance(self, L, x, b):
        raise NotImplementedError

    def _adaptRHSTolerance(self, L, x, b):
        raise NotImplementedError

    def _adaptMatrixTolerance(self, L, x, b):
        raise NotImplementedError

    def _adaptInitialTolerance(self, L, x, b):
        raise NotImplementedError

    def _adaptSolutionTolerance(self, L, x, b):
        raise NotImplementedError

    def _adaptPreconditionedTolerance(self, L, x, b):
        raise NotImplementedError

    def _adaptNaturalTolerance(self, L, x, b):
        raise NotImplementedError

    def _adaptTolerance(self, L, x, b):
        adapt = {
            "default": self._adaptDefaultTolerance,
            "unscaled": self._adaptUnscaledTolerance,
            "RHS": self._adaptRHSTolerance,
            "matrix": self._adaptMatrixTolerance,
            "initial": self._adaptInitialTolerance,
            "solution": self._adaptSolutionTolerance,
            "preconditioned": self._adaptPreconditionedTolerance,
            "natural": self._adaptNaturalTolerance
        }

        tolerance_scale, suite_criterion = adapt[self.criterion](L, x, b)

        return tolerance_scale, suite_criterion

    @staticmethod
    def scale_tolerance(tol, scale):
        if tol is not None:
            tol *= scale
        return tol

    def _applyUnderRelaxation(self, underRelaxation=None):
        if underRelaxation is not None:
            self.matrix.putDiagonal(numerix.asarray(self.matrix.takeDiagonal()) / underRelaxation)
            self.RHSvector += (1 - underRelaxation) * self.matrix.takeDiagonal() * numerix.array(self.var).flatten()

    def _calcResidualVector(self, residualFn=None):
        if residualFn is not None:
            return residualFn(self.var, self.matrix, self.RHSvector)
        else:
            Lx = self.matrix * numerix.array(self.var).flatten()

            return Lx - self.RHSvector

    def _calcResidual(self, residualFn=None):
        if residualFn is not None:
            return residualFn(self.var, self.matrix, self.RHSvector)
        else:
            return numerix.L2norm(self._calcResidualVector())

    def _calcRHSNorm(self):
        return numerix.L2norm(self.RHSvector)

    def __repr__(self):
        return '%s(tolerance=%g, iterations=%g)' \
            % (self.__class__.__name__, self.tolerance, self.iterations)

    def _canSolveAsymmetric(self):
        return True

    def __enter__(self):
        # The __enter__() and __exit__() methods allow
        # solver objects to be used within a context manager.
        # This enables one to write:
        #
        #     with ExampleSolver() as solver:
        #         eq.solve(var=phi, solver=solver)
        #
        # For subclasses that do not define __enter__()
        # and __exit__() methods, the above is equivalent to:
        #
        #     eq.solve(var=phi, solver=ExampleSolver())

        return self

    def __exit__(self, exc_type, exc_value, traceback):
        pass


    def __del__(self):
        pass

    def _test(self):
        """
        >>> import fipy as fp

        Consider a steady-state 1D diffusion problem with a
        position-dependent diffusivity and Dirichlet boundary conditions:

        .. math::

           \begin{aligned}
           \frac{\partial}{\partial x}\left[
               \left(1 + x\right)
               \frac{\partial \phi}{\partial x}
           \right] &= 0
           \\
           \left.\phi\right\rvert_{x=0} &= \phi_L
           \\
           \left.\phi\right\rvert_{x=1} &= \phi_R
           \end{aligned}

        with the analytical solution

        .. math::

           \phi = \frac{\phi_R - \phi_L}{\ln 2} \ln\left(1 + x\right) + \phi_L

        >>> N = 100
        >>> mesh = fp.Grid1D(nx=N, Lx=1)
        >>> phi = fp.CellVariable(mesh=mesh, name=r"$\phi")
        >>> phiL = 1000.
        >>> phiR = 2000.
        >>> phi_analytical = ((((phiR - phiL)/fp.numerix.log(2.))
        ...                    * fp.numerix.log(1 + mesh.x))
        ...                   + phiL)
        >>> phi_analytical.name = r"$\phi_{analytical}$"

        >>> fp.numerix.random.seed(12345)
        >>> phi_initial = phi_analytical + fp.GaussianNoiseVariable(mesh=mesh, variance=1e-3)
        >>> phi.value = phi_initial
        >>> phi.constrain(phiL, where=mesh.facesLeft)
        >>> phi.constrain(phiR, where=mesh.facesRight)
        >>> D = fp.FaceVariable(mesh=mesh, value=1 + mesh.faceCenters[0])
        >>> eq = fp.DiffusionTerm(coeff=D) == 0

        The norm of the matrix will be :math:`\mathcal{O}(8 N)`.
        The norm of the right-hand side will be :math:`math{O}(8000 N}`.
        We choose the initial condition such that the order of the initial
        residual is :math:`\mathcal{O}(400 / N}`.

        >>> # Solver = fp.LinearPCGSolver
        >>> # Solver = fp.LinearGMRESSolver
        >>> # Solver = fp.LinearGMRESSolver
        >>> Solver = fp.LinearCGSSolver

        >>> solver = Solver(precon=None)
        >>> solver = eq._prepareLinearSystem(var=phi, solver=solver, boundaryConditions=(), dt=1.)
        >>> Lnorm, bnorm, rnorm = solver._norms
        >>> enorm = fp.numerix.L2norm(phi - phi_analytical) / fp.numerix.L2norm(phi_analytical)
        >>> print("|L| = {Lnorm}".format(**locals()))
        >>> print("|b| = {bnorm}".format(**locals()))
        >>> print("|r| = {rnorm}".format(**locals()))
        >>> print("|e| = {enorm}".format(**locals()))

        >>> criteria = [
        ...     ("unscaled", 1.),
        ...     ("RHS", bnorm),
        ...     ("matrix", Lnorm),
        ...     ("initial", rnorm)
        ... ]
        >>> # criteria += ["solution"]  doctest: +TRILINOS_SOLVER
        >>> criteria += [
        ...     ("preconditioned", bnorm),
        ...     ("natural", bnorm)
        ... ] # doctest: +PETSC_SOLVER
        >>> for (criterion, target) in criteria:
        ...     phi.setValue(phi_initial)
        ...     with Solver(criterion=criterion, precon=None, tolerance=1e-5) as s:
        ...         res = eq.sweep(var=phi, solver=s)
        ...         # print(s.convergence)
        ...         print(",".join([s.convergence.suite,
        ...                         criterion,
        ...                         s.convergence.status_name,
        ...                         str(s.convergence.residual),
        ...                         str(s.convergence.residual / (s.tolerance * target)),
        ...                         str(s.convergence.iterations),
        ...                         str(fp.numerix.L2norm(phi - phi_analytical) / fp.numerix.L2norm(phi_analytical))
        ...                        ]))
        """
        pass

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
