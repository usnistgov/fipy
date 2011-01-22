from smoothedAggregationSolver import SmoothedAggregationSolver
from generalSolver import GeneralSolver
from fipy.solvers.pysparse.linearPCGSolver import LinearPCGSolver

"""
try:
    from pyamg import solveit
    DefaultSolver = GeneralSolver
    DefaultAsymmetricSolver = GeneralSolver
except ImportError:
"""
DefaultSolver = SmoothedAggregationSolver
DefaultAsymmetricSolver = SmoothedAggregationSolver
