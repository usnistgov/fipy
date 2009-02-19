from fipy.solvers.solver import Solver

from pykrylov.cg import CG
from pykrylov.cgs import CGS
from fipy.tools.pysparseMatrix import _PysparseMatrix
from fipy.tools import numerix


class PykrylovSolver(Solver):
    def _solve(self, L, x, b):

        diag = L.takeDiagonal()

        jacobi =  self._getMatrixClass()(L.matrix.shape[0])
        jacobi.putDiagonal(1. / diag)
        L = jacobi * L
        b = b / diag

        solver = self.solverClass(lambda u: L * u, reltol=self.tolerance, abstol=self.tolerance, matvec_max=self.iterations)

        solver.solve(b, guess=x)
        x[:] = solver.bestSolution
##        print 'solver.nMatvec',solver.nMatvec
##        print 'solver.residNorm',solver.residNorm


    def _getMatrixClass(self):
        return _PysparseMatrix
        
class LinearCGSSolver(PykrylovSolver):
    solverClass = CGS

    def _canSolveAssymetric(self):
        return True
    
class LinearPCGSolver(PykrylovSolver):
    solverClass = CG

    def _canSolveAssymetric(self):
        return False

LinearPCGSolver = LinearCGSSolver
DefaultSolver = LinearCGSSolver
from fipy.solvers.pysparse.linearLUSolver import LinearLUSolver
##LinearLUSolver = LinearCGSSolver

