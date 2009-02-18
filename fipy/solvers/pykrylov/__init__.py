from fipy.solvers.solver import Solver

from pykrylov.cg import CG
from pykrylov.cgs import CGS
from fipy.tools.pysparseMatrix import _PysparseMatrix
from fipy.tools import numerix


class PykrylovSolver(Solver):
    def _solve(self, L, x, b):
        diagL = numerix.absolute(L.takeDiagonal())
        solver = self.solverClass(lambda u: L * u, reltol=self.tolerance, matvec_max=self.iterations)##, precon=lambda u: u / diagL)
        solver.solve(b, guess=x)
        x[:] = solver.bestSolution

    def _getMatrixClass(self):
        return _PysparseMatrix
        
class LinearCGSSolver(PykrylovSolver):
    solverClass = CGS

class LinearPCGSolver(PykrylovSolver):
    solverClass = CG

DefaultSolver = LinearPCGSolver
LinearGMRESSolver = LinearCGSSolver
LinearLUSolver = LinearCGSSolver
    
