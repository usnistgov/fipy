#!/usr/bin/env python

##
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "solver.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this software is not subject to copyright
 # protection and is in the public domain.  FiPy is an experimental
 # system.  NIST assumes no responsibility whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 #
 # This software can be redistributed and/or modified freely
 # provided that any derivative works bear some notice that they are
 # derived from it, and any modified versions bear some notice that
 # they have been modified.
 # ========================================================================
 #
 # ###################################################################
 ##

"""
The iterative solvers may output warnings if the solution is considered
unsatisfactory. If you are not interested in these warnings, you can invoke
python with a warning filter such as::

    $ python -Wignore::fipy.SolverConvergenceWarning myscript.py

If you are extremely concerned about your preconditioner for some reason, you
can abort whenever it has problems with::

    $ python -Werror::fipy.PreconditionerWarning myscript.py

"""
__docformat__ = 'restructuredtext'

from fipy.tools import numerix

__all__ = ["SolverConvergenceWarning", "MaximumIterationWarning",
           "PreconditionerWarning", "IllConditionedPreconditionerWarning",
           "PreconditionerNotPositiveDefiniteWarning", "MatrixIllConditionedWarning",
           "StagnatedSolverWarning", "ScalarQuantityOutOfRangeWarning", "Solver"]

class SolverConvergenceWarning(Warning):
    def __init__(self, solver, iter, relres):
        self.solver = solver
        self.iter = iter
        self.relres = relres

    def __str__(self):
        return "%s failed. Iterations: %g. Relative error: %g" % (str(self.solver), self.iter, self.relres)

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

class Solver(object):
    """
    The base `LinearXSolver` class.

    .. attention:: This class is abstract. Always create one of its subclasses.
    """

    def __init__(self, tolerance=1e-10, iterations=1000, precon=None):
        """
        Create a `Solver` object.

        :Parameters:
          - `tolerance`: The required error tolerance.
          - `iterations`: The maximum number of iterative steps to perform.
          - `precon`: Preconditioner to use. This parameter is only available for Trilinos solvers.

        """
        if self.__class__ is Solver:
            raise NotImplementedError("can't instantiate abstract base class")

        self.tolerance = tolerance
        self.iterations = iterations

        self.preconditioner = precon

    def _storeMatrix(self, var, matrix, RHSvector):
        self.var = var
        self.matrix = matrix
        self.RHSvector = RHSvector

    def _solve(self):
        raise NotImplementedError

    def _solve_(self, L, x, b):
        raise NotImplementedError

    def _applyUnderRelaxation(self, underRelaxation=None):
        if underRelaxation is not None:
            self.matrix.putDiagonal(self.matrix.takeDiagonal() / underRelaxation)
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

    _warningList = (ScalarQuantityOutOfRangeWarning,
                    StagnatedSolverWarning,
                    MatrixIllConditionedWarning,
                    PreconditionerNotPositiveDefiniteWarning,
                    IllConditionedPreconditionerWarning,
                    MaximumIterationWarning)

    def _raiseWarning(self, info, iter, relres):
        # info is negative, so we list in reverse order so that
        # info can be used as an index from the end

        if info < 0:
            # is stacklevel=5 always what's needed to get to the user's scope?
            import warnings
            warnings.warn(self._warningList[info](self, iter, relres), stacklevel=5)

    def __repr__(self):
        return '%s(tolerance=%g, iterations=%g)' \
            % (self.__class__.__name__, self.tolerance, self.iterations)

    def _canSolveAsymmetric(self):
        return True
