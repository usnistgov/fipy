import os
import sys

from solver import SolverConvergenceWarning, \
     PreconditionerWarning, \
     ScalarQuantityOutOfRangeWarning, \
     StagnatedSolverWarning, \
     MatrixIllConditionedWarning, \
     PreconditionerNotPositiveDefiniteWarning, \
     IllConditionedPreconditionerWarning, \
     MaximumIterationWarning

args = [s.lower() for s in sys.argv[1:]]

# any command-line specified solver takes precedence over environment variables
if '--trilinos' in args:
    solver = "trilinos"
elif '--pysparse' in args:
    solver = "pysparse"
elif '--no-pysparse' in args:
    solver = "no-pysparse"
elif os.environ.has_key('FIPY_SOLVERS'):
    solver = os.environ['FIPY_SOLVERS'].lower()
else:
    solver = None
    
if solver == "pysparse":
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
elif solver == "no-pysparse":
    from fipy.solvers.trilinos import *
    from fipy.matrices.trilinosMatrix import _TrilinosMeshMatrix
    _MeshMatrix =  _TrilinosMeshMatrix 
elif solver is None:
    # If no argument or environment variable, try importing them and seeing
    # what works
    try: 
        from fipy.tools import parallel
        if parallel.Nproc > 1:
            raise Exception("PySparse solvers cannot be used with multiple processors")

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
                from fipy.matrices.trilinosMatrix import _TrilinosMeshMatrix
                _MeshMatrix =  _TrilinosMeshMatrix
        except:
            raise ImportError, "Could not import any solver package. If you are using Trilinos, make sure you have all of the necessary Trilinos packages installed - Epetra, EpetraExt, AztecOO, Amesos, ML, and IFPACK." 
else:
    raise ImportError, 'Unknown solver package %s' % solver

