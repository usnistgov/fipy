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
import os
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

    #: Default tolerance for linear solves unless `criterion="legacy"`
    DEFAULT_TOLERANCE = 1e-5

    #: Default tolerance for linear solves if `criterion="legacy"`
    LEGACY_TOLERANCE = 1e-10

    #: Default maximum number of iterative steps to perform
    DEFAULT_ITERATIONS = 1000

    #: Default preconditioner to apply to the matrix
    DEFAULT_PRECONDITIONER = None

    def __init__(self, tolerance="default", criterion="default",
                 iterations="default", precon="default"):
        """
        Create a `Solver` object.

        Parameters
        ----------
        tolerance : float
            Required residual tolerance.
        criterion : {'default', 'initial', 'unscaled', 'RHS', 'matrix', 'solution', 'preconditioned', 'natural', 'legacy'}, optional
            Interpretation of ``tolerance``.
            See :ref:`CONVERGENCE` for more information.
        iterations : int
            Maximum number of iterative steps to perform.
        precon : ~fipy.solvers.preconditioner.Preconditioner
            Preconditioner to use.  Not all solver suites support
            preconditioners.
        """
        if self.__class__ is Solver:
            raise NotImplementedError("can't instantiate abstract base class")

        self.criterion = self.value_or_default(criterion,
                                               "default",
                                               "FIPY_DEFAULT_CRITERION")
        self.tolerance = self.value_or_default(tolerance,
                                               self.default_tolerance)
        self.iterations = self.value_or_default(iterations,
                                                self.DEFAULT_ITERATIONS)
        self.preconditioner = self.value_or_default(precon,
                                                    self.default_preconditioner)

        self._log = logging.getLogger(self.__class__.__module__
                                      + "." + self.__class__.__name__)

    def value_or_default(self, value, default, envvar=None):
        if value == "default":
            if envvar is not None:
                value = os.environ.get(envvar, default)
            else:
                value = default
        return value

    @property
    def default_tolerance(self):
        """Default tolerance for linear solve
        """
        if self.criterion == "legacy":
            return self.LEGACY_TOLERANCE
        else:
            return self.DEFAULT_TOLERANCE

    @property
    def default_preconditioner(self):
        if self.DEFAULT_PRECONDITIONER is not None:
            # instantiate DEFAULT_PRECONDITIONER class
            return self.DEFAULT_PRECONDITIONER()
        else:
            return None

    def _storeMatrix(self, var, matrix, RHSvector):
        self.var = var
        self.matrix = matrix
        self.RHSvector = RHSvector

    def _scatterGhosts(self, x):
        """Distribute ghost values (if any) across processes
        """
        return x

    def _cleanup(self):
        pass

    def _solve(self):
        """Solve system of equations posed for FiPy

        Common method invoked by :class:`~fipy.terms.term.Term`, which then
        calls solver-suite-specific :math:`~fipy.solvers.solver._solve_`
        methods.
        """
        L, x, b = self._Lxb

        if self._log.isEnabledFor(logging.DEBUG):
            s = "Lnorm: {Lnorm}, bnorm: {bnorm}, rnorm: {rnorm}"
            self._log.debug(s.format(Lnorm=self._matrixNorm(L, x, b),
                                     bnorm=self._rhsNorm(L, x, b),
                                     rnorm=self._residualNorm(L, x, b)))

        x = self._solve_(L, x, b)

        x = self._scatterGhosts(x)

        factor = self.var.unit.factor
        if factor != 1:
            x /= self.var.unit.factor

        self.var.value = x.reshape(self.var.shape)

        self._cleanup()

    def _solve_(self, L, x, b):
        """Solve system of equations posed for solver suite

        Parameters
        ----------
        L : ~fipy.matrices.sparseMatrix._SparseMatrix
            Sparse matrix object
        x : array_like
            Solution variable in form suitable for solver
        b : array_like
            Right-hand side vector in form suitable for solver

        Returns
        -------
        ndarray
            Solution vector
        """
        raise NotImplementedError

    def _setConvergence(self, suite, code, iterations, residual, actual_code=None, **kwargs):
        cls = ConvergenceBase.code_registry[(suite, code)]
        self.convergence = cls(solver=self,
                               iterations=iterations,
                               residual=residual,
                               criterion=self.criterion,
                               actual_code=actual_code,
                               **kwargs)

    def _legacyNorm(self, L, x, b):
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
        """Matrix, solution vector, and right-hand side vector

        Returns
        -------
        L : matrix
            Sparse matrix in form suitable for solver
        x : ndarray
            Solution variable in form suitable for solver
        b : ndarray
            Right-hand side vector in form suitable for solver
        """
        if self.var.mesh.communicator.Nproc > 1:
            raise Exception(str(type(self)) + " cannot be used with multiple processors")

        L = self.matrix.matrix
        x = self.var.numericValue.ravel()
        b = numerix.asarray(self.RHSvector)

        if ((self.matrix == 0)
            or (L.shape[0] != L.shape[1])
            or (L.shape[0] != len(x))):

            from fipy.terms import SolutionVariableNumberError

            raise SolutionVariableNumberError

        return (L, x, b)

    def _adaptLegacyTolerance(self, L, x, b):
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
            "legacy": self._adaptLegacyTolerance,
            "unscaled": self._adaptUnscaledTolerance,
            "RHS": self._adaptRHSTolerance,
            "matrix": self._adaptMatrixTolerance,
            "initial": self._adaptInitialTolerance,
            "solution": self._adaptSolutionTolerance,
            "preconditioned": self._adaptPreconditionedTolerance,
            "natural": self._adaptNaturalTolerance,
            "default": self._adaptRHSTolerance
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
        r"""
        >>> import fipy as fp

        For sufficiently constrained circumstances, all solver suites
        should do the same thing.  The following problem setup is designed
        to ensure that all interpret solver criteria correctly and achieve
        the "same" tolerance in the same number of iterations.

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
        >>> phi_analytical.name = r"$\phi_\mathrm{analytical}$"

        >>> fp.numerix.random.seed(12345)
        >>> variance = 1e-3
        >>> phi_initial = phi_analytical + fp.GaussianNoiseVariable(mesh=mesh, variance=variance)
        >>> phi.value = phi_initial
        >>> phi.constrain(phiL, where=mesh.facesLeft)
        >>> phi.constrain(phiR, where=mesh.facesRight)
        >>> D = fp.FaceVariable(mesh=mesh, value=1 + mesh.faceCenters[0])
        >>> eq = fp.DiffusionTerm(coeff=D) == 0

        For reproducibility between suites, we select a solver with
        predictable characteristics (that counts out GMRES) and no
        preconditioning.

        >>> Solver = fp.LinearCGSSolver
        >>> solver = Solver(precon=None)

        >>> solver = eq._prepareLinearSystem(var=phi,
        ...                                  solver=solver,
        ...                                  boundaryConditions=(),
        ...                                  dt=1.)
        >>> L, x, b = solver._Lxb

        The problem parameters were chosen to give good separation between the
        different convergence norms.

        The norm of the matrix is the infinity norm

        .. math::

           \left\| L_{ij}\right\|_\infty &= \max_i \sum_j \left| A_ij \right|
           \\
           &= \max_i \left[
               \left| -N(1 + x_i) \right|
               + \left| 2N(1 + x_i) \right|
               + \left| -N(1 + x_i) \right|
           \right]
           \\
           &= \max_i 4N(1 + x_i)
           &= \mathcal{O}(8 N)

        >>> Lnorm = solver._matrixNorm(L, x, b)
        >>> print(numerix.allclose(Lnorm, 8 * N, rtol=0.1))
        True

        The right-hand-side vector is zero except at the boundaries,
        where the contribution is

        .. math::

           \frac{(1 + x) \phi_{BC} A_f}{d_{AP}} &= (1 + x) \phi_{BC} 2 N
           \\
           &= 2 N \phi_L = 2000 N\qquad\text{at $x = 0$}
           \\
           &= 4 N \phi_R = 8000 N\qquad\text{at $x = 1$}

        Thus the :math:`L_2` norm of the right-hand-side vector is
        :math:`\left\| b \right\|_2 = \math{O}(8000 N}`.

        >>> bnorm = solver._rhsNorm(L, x, b)
        >>> print(numerix.allclose(bnorm, 8000 * N, rtol=0.1))
        True

        We choose the initial condition such that the initial residual will
        be small.

        .. math::

           \phi_0 &= \phi_\text{analytical} + \mathcal{O}(\sigma)
           \\
           r = L \phi_0 - b
           &= L \phi_\text{analytical} - b + L \mathcal{O}(\sigma)
           \\
           &= L \mathcal{O}(\sigma)
           \\
           \left\| r \right\|_2 &= \left\| L \mathcal{O}(\sigma) \right\|_2
           \\
           &= \sqrt{\sum_{0 \le i < N} \left[
               N(1 + x_i) \mathcal{O}(\sigma)
               + 2N(1 + x_i) \mathcal{O}(\sigma)
               + N(1 + x_i) \mathcal{O}(\sigma)
           \right]^2}
           \\
           &= 4 N \mathcal{O}(\sigma) \sqrt{\sum_{0 \le i < N} (1 + x_i)^2}
           \\
           &= \text{probably $\sqrt{\pi}$ or something}
           \\
           &= \mathcal{O}(4 N \sqrt{N} \sigma)

        >>> rnorm = solver._residualNorm(L, x, b)
        >>> print(numerix.allclose(rnorm, 4 * N * numerix.sqrt(N * variance),
        ...                        rtol=0.1))
        True

        Calculate the error of the initial condition (probably could be
        estimated via truncation error blah blah blah).

        >>> enorm = fp.numerix.L2norm(phi - phi_analytical) / fp.numerix.L2norm(phi_analytical)

        >>> from fipy.solvers.convergence import Convergence

        Check that:
        - the solution is converged,
        - the solver reaches the desired residual for the
          criterion, without overshooting too much.  Most get close, but
          "unscaled" overshoots a lot for most suites.
        - the iteration count is as expected
        - the error has been reduced from the initial guess

        SciPy 1.12 translated all of their solvers from FORTRAN to Python
        and the iteration counts went up for some reason.

        >>> criteria = [
        ...     ("unscaled", 1., 0.003, 118),
        ...     ("RHS", bnorm, 0.6, 2),
        ...     ("matrix", Lnorm, 0.6, 60),
        ...     ("initial", rnorm, 0.6, 114)
        ... ] # doctest: +SCIPY_PYTHON_SOLVER
        >>> criteria = [
        ...     ("unscaled", 1., 0.003, 114),
        ...     ("RHS", bnorm, 0.6, 2),
        ...     ("matrix", Lnorm, 0.6, 58),
        ...     ("initial", rnorm, 0.6, 110)
        ... ] # doctest: +NOT_SCIPY_PYTHON_SOLVER

        >>> # criteria += ["solution"]  doctest: +TRILINOS_SOLVER
        >>> criteria += [
        ...     ("preconditioned", bnorm, 0.6, 2),
        ...     ("natural", bnorm, 0.6, 6)
        ... ] # doctest: +PETSC_SOLVER
        >>> satisfied = []
        >>> for (criterion, target, lower_bound, iterations) in criteria:
        ...     phi.setValue(phi_initial)
        ...     with Solver(criterion=criterion, precon=None) as s:
        ...         res = eq.sweep(var=phi, solver=s)
        ...         error = (fp.numerix.L2norm(phi - phi_analytical)
        ...                  / fp.numerix.L2norm(phi_analytical))
        ...         checks = [isinstance(s.convergence, Convergence),
        ...                   (lower_bound
        ...                    < (s.convergence.residual
        ...                       / (s.tolerance * target))
        ...                    < 1.0),
        ...                   numerix.allclose(s.convergence.iterations,
        ...                                    iterations,
        ...                                    atol=1),
        ...                   error < enorm]
        ...         print(criterion, s.convergence, target, lower_bound, s.convergence.residual / (s.tolerance * target), iterations, s.convergence.iterations, error, enorm)
        ...         satisfied.append(all(checks))
        >>> print(all(satisfied))
        True
        >>> print(satisfied)
        """
        pass

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
