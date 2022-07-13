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

import json
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
        criterion : {'default', 'initial', 'unscaled', 'RHS', 'matrix', 'solution', 'preconditioned', 'natural'}
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

        info = self.convergence.info.copy()
        info["solver"] = str(info["solver"])
        self._log.debug(json.dumps(info))

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

        tolerance_factor, suite_criterion = adapt[self.criterion](L, x, b)

        return tolerance_factor, suite_criterion

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

    def _test(self):
        """
        >>> import fipy as fp

        >>> mesh = fp.Grid1D(nx=3)
        >>> var = fp.CellVariable(mesh=mesh)
        >>> var.constrain(1., where=mesh.facesLeft)
        >>> var.constrain(2., where=mesh.facesRight)
        >>> D = fp.FaceVariable(mesh=mesh, value=mesh.faceCenters[0])
        >>> eq = fp.TransientTerm() == fp.DiffusionTerm(coeff=D)

        >>> var.setValue(mesh.x)
        >>> # with eq.getDefaultSolver(criterion="default", precon=None) as s:
        >>> # with fp.LinearPCGSolver(criterion="default") as solver:
        >>> # with fp.LinearGMRESSolver(criterion="default") as solver:
        >>> with eq.getDefaultSolver(criterion="default") as s:
        ...     eq.solve(var=var, dt=1., solver=s)
        ...     print(s.convergence.criterion, s.convergence.residual, var)
        default 5.325944231618159e-16 [ 1.06363636  1.62727273  1.97272727]

        >>> var.setValue(mesh.x)
        >>> with eq.getDefaultSolver(criterion="unscaled") as s:
        ...     eq.solve(var=var, dt=1., solver=s)
        ...     print(s.convergence.criterion, s.convergence.residual, var)
        unscaled 8.479955780533376e-15 [ 1.06363636  1.62727273  1.97272727]

        >>> var.setValue(mesh.x)
        >>> with eq.getDefaultSolver(criterion="RHS") as s:
        ...     eq.solve(var=var, dt=1., solver=s)
        ...     print(s.convergence.criterion, s.convergence.residual, var)
        RHS 5.813782806744329e-16 [ 1.06363636  1.62727273  1.97272727]

        >>> var.setValue(mesh.x)
        >>> with eq.getDefaultSolver(criterion="matrix") as s:
        ...     eq.solve(var=var, dt=1., solver=s)
        ...     print(s.convergence.criterion, s.convergence.residual, var)
        matrix 7.709050709575797e-16 [ 1.06363636  1.62727273  1.97272727]

        >>> var.setValue(mesh.x)
        >>> with eq.getDefaultSolver(criterion="initial") as s:
        ...     eq.solve(var=var, dt=1., solver=s)
        ...     print(s.convergence.criterion, s.convergence.residual, var)
        initial 1.6319682508690227e-15 [ 1.06363636  1.62727273  1.97272727]

        >>> # var.setValue(mesh.x)
        >>> # with eq.getDefaultSolver(criterion="solution") as s: # doctest: +TRILINOS_SOLVER
        ... #     eq.solve(var=var, dt=1., solver=s)
        ... #     print(s.convergence.criterion, s.convergence.residual, var)

        >>> var.setValue(mesh.x)
        >>> with eq.getDefaultSolver(criterion="preconditioned") as s: # doctest: +PETSC_SOLVER
        ...     eq.solve(var=var, dt=1., solver=s)
        ...     print(s.convergence.criterion, s.convergence.residual, var)
        preconditioned 5.325944231618159e-16 [ 1.06363636  1.62727273  1.97272727]

        >>> var.setValue(mesh.x)
        >>> with eq.getDefaultSolver(criterion="natural") as s: # doctest: +PETSC_SOLVER
        ...     eq.solve(var=var, dt=1., solver=s)
        ...     print(s.convergence.criterion, s.convergence.residual, var)
        """
        pass

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
