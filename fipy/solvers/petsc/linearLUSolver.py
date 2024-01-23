from __future__ import division
from builtins import range
from past.utils import old_div
__docformat__ = 'restructuredtext'

from petsc4py import PETSc

from fipy.solvers.petsc.petscSolver import PETScSolver
from fipy.tools.timer import Timer

__all__ = ["LinearLUSolver"]

class LinearLUSolver(PETScSolver):

    """
    The `LinearLUSolver` is an interface to the LU preconditioner in PETSc.
    A direct solve is performed.

    """

    def __init__(self, tolerance=1e-10, iterations=10, precon="lu"):
        """
        :Parameters:
          - `tolerance`: The required error tolerance.
          - `iterations`: The maximum number of iterative steps to perform.
          - `precon`: *Ignored*.

        """
        PETScSolver.__init__(self, tolerance=tolerance,
                             iterations=iterations, precon="lu")

    def _solve_(self, L, x, b):
        ksp = PETSc.KSP()
        ksp.create(PETSc.COMM_WORLD)
        ksp.setType("preonly")
        ksp.getPC().setType(self.preconditioner)
        # TODO: SuperLU invoked with PCFactorSetMatSolverType(pc, MATSOLVERSUPERLU)
        #       see: http://www.mcs.anl.gov/petsc/petsc-dev/src/ksp/ksp/examples/tutorials/ex52.c.html
        # PETSc.PC().setFactorSolverType("superlu")

        L.assemble()
        ksp.setOperators(L)
        ksp.setFromOptions()

        self._log.debug("BEGIN solve")

        with Timer() as t:
            for iteration in range(self.iterations):
                errorVector = L * x - b
                tol = errorVector.norm()

                if iteration == 0:
                    tol0 = tol

                if tol <= self.tolerance * tol0:
                    break

                xError = x.copy()

                ksp.solve(errorVector, xError)
                x -= xError

        self._log.debug("END solve - {} ns".format(t.elapsed))

        self._log.debug('solver: %s', ksp.type)
        self._log.debug('precon: %s', ksp.getPC().type)
        self._log.debug('iterations: %d / %d', iteration+1, self.iterations)
        self._log.debug('residual: %s', errorVector.norm(1))
