import sys
from solver import SolverConvergenceWarning, \
     PreconditionerWarning, \
     ScalarQuantityOutOfRangeWarning, \
     StagnatedSolverWarning, \
     MatrixIllConditionedWarning, \
     PreconditionerNotPositiveDefiniteWarning, \
     IllConditionedPreconditionerWarning, \
     MaximumIterationWarning

# First check for command-line arguments
if '--Trilinos' in sys.argv[1:]:
    from fipy.solvers.trilinos import *
elif '--Pysparse' in sys.argv[1:]:
    from fipy.solvers.pysparse import *
else:
    import os
    # Next, check for an environment variable telling us which solver to use
    if os.environ.has_key('FIPY_SOLVERS'):
        if os.environ['FIPY_SOLVERS'].lower() == 'pysparse':
            from fipy.solvers.pysparse import *
        elif os.environ['FIPY_SOLVERS'].lower() == 'trilinos':
            from fipy.solvers.trilinos import *
        else:
            raise ImportError, 'Unknown solver package %s' % os.environ['FIPY_SOLVERS']
    else:
        # If no argument or environment variable, try importing them and seeing
        # what works
        if not foundSolvers:
            try: 
                from fipy.solvers.pysparse import *
            except:
                try:
                    from fipy.solvers.trilinos import *
                except:
                    raise ImportError, "Could not import any solver package. If you are using Trilinos, make sure you have all of the necessary Trilinos packages installed - Epetra, EpetraExt, AztecOO, Amesos, ML, and IFPACK." 
