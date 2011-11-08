from fipy.tools.parser import _parseSolver
from fipy.tools  import parallel as _parallel

from solver import SolverConvergenceWarning, \
     PreconditionerWarning, \
     ScalarQuantityOutOfRangeWarning, \
     StagnatedSolverWarning, \
     MatrixIllConditionedWarning, \
     PreconditionerNotPositiveDefiniteWarning, \
     IllConditionedPreconditionerWarning, \
     MaximumIterationWarning

solver = _parseSolver()

def _envSolver(solver):
    import os
    if solver is None and os.environ.has_key('FIPY_SOLVERS'):
        solver = os.environ['FIPY_SOLVERS'].lower()
    return solver
    
solver = _envSolver(solver)

if solver == "pysparse":
    if _parallel.Nproc > 1:
        raise  Exception('pysparse solvers do not run in parallel')
    from fipy.solvers.pysparse import *
    from fipy.matrices.pysparseMatrix import _PysparseMeshMatrix
    _MeshMatrix =  _PysparseMeshMatrix

elif solver == "trilinos":
    from fipy.solvers.trilinos import *

    try:
        from fipy.matrices.pysparseMatrix import _PysparseMeshMatrix
        _MeshMatrix =  _PysparseMeshMatrix
    except ImportError:
        from fipy.matrices.trilinosMatrix import _TrilinosMeshMatrix
        _MeshMatrix =  _TrilinosMeshMatrix

elif solver == "scipy":
    if _parallel.Nproc > 1:
        raise  Exception('scipy solvers do not run in parallel')
    from fipy.solvers.scipy import *
    from fipy.matrices.scipyMatrix import _ScipyMeshMatrix
    _MeshMatrix = _ScipyMeshMatrix
    
elif solver == "pyamg":
    if _parallel.Nproc > 1:
        raise  Exception('pyamg solvers do not run in parallel')
    from fipy.solvers.pyAMG import *
    from fipy.matrices.scipyMatrix import _ScipyMeshMatrix
    _MeshMatrix = _ScipyMeshMatrix
    
elif solver == "no-pysparse":
    from fipy.solvers.trilinos import *
    from fipy.matrices.trilinosMatrix import _TrilinosMeshMatrix
    _MeshMatrix =  _TrilinosMeshMatrix 

elif solver is None:
    # If no argument or environment variable, try importing them and seeing
    # what works

    
   
    try:
        if _parallel.Nproc > 1:
            raise  Exception('pysparse solvers do not run in parallel')
        from fipy.solvers.pysparse import *
        solver = "pysparse"
        from fipy.matrices.pysparseMatrix import _PysparseMeshMatrix
        _MeshMatrix =  _PysparseMeshMatrix
    except:
        try:
            from fipy.solvers.trilinos import *
            solver = "trilinos"
            try:
                from fipy.matrices.pysparseMatrix import _PysparseMeshMatrix
                _MeshMatrix =  _PysparseMeshMatrix
            except ImportError:
                solver = "no-pysparse"
                from fipy.matrices.trilinosMatrix import _TrilinosMeshMatrix
                _MeshMatrix =  _TrilinosMeshMatrix
        except:
            try:
                if _parallel.Nproc > 1:
                    raise  Exception('pyamg solvers do not run in parallel')
                from fipy.solvers.pyAMG import *
                solver = "pyamg"
                from fipy.matrices.scipyMatrix import _ScipyMeshMatrix
                _MeshMatrix = _ScipyMeshMatrix
            except:
                try:
                    if _parallel.Nproc > 1:
                        raise  Exception('scipy solvers do not run in parallel')
                    from fipy.solvers.scipy import *
                    solver = "scipy"
                    from fipy.matrices.scipyMatrix import _ScipyMeshMatrix
                    _MeshMatrix = _ScipyMeshMatrix
                except:
                    raise ImportError, "Could not import any solver package. If you are using Trilinos, make sure you have all of the necessary Trilinos packages installed - Epetra, EpetraExt, AztecOO, Amesos, ML, and IFPACK." 
else:
    raise ImportError, 'Unknown solver package %s' % solver

