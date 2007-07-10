import sys
from solver import SolverConvergenceWarning, \
     PreconditionerWarning, \
     ScalarQuantityOutOfRangeWarning, \
     StagnatedSolverWarning, \
     MatrixIllConditionedWarning, \
     PreconditionerNotPositiveDefiniteWarning, \
     IllConditionedPreconditionerWarning, \
     MaximumIterationWarning

if '--Trilinos' in sys.argv[1:]:
    from fipy.solvers.trilinos import *
elif '--Pysparse' in sys.argv[1:]:
    from fipy.solvers.pysparse import *
elif '--Scipy' in sys.argv[1:]:
    from fipy.solvers.scipy import *
else:
    import os
    if os.environ.has_key('FIPY_SOLVERS'):
        if os.environ['FIPY_SOLVERS'] == 'Pysparse':
            from fipy.solvers.pysparse import *
        elif os.environ['FIPY_SOLVERS'] == 'Trilinos':
            from fipy.solvers.trilinos import *
        elif os.environ['FIPY_SOLVERS'] == 'Scipy':
            from fipy.solvers.scipy import *
        else:
            raise ImportError, 'Unknown solver package %s' % os.environ['FIPY_SOLVERS']
    else:
        foundSolvers  = False
        if not foundSolvers:
            try: 
                from fipy.solvers.pysparse import *
                foundSolvers = True
            except:
                pass

        if not foundSolvers:
            try:
                from fipy.solvers.trilinos import *
                foundSolvers = True
            except:
                pass
        
        if not foundSolvers:
            try: 
                from fipy.solvers.scipy import *
                foundSolvers = True
            except:
                pass

        if not foundSolvers:
            raise ImportError, "Could not import any solver package. If you are using Trilinos, make sure you have all of the necessary Trilinos packages installed - Epetra, EpetraExt, AztecOO, Amesos, ML, and IFPACK." 

if(solverSuite() != 'Trilinos'):
    def mainProcessor():
        return True
else:
    from PyTrilinos import Epetra
    if(Epetra.PyComm().MyPID() == 0):
        def mainProcessor():
            return True
    else:
        def mainProcessor():
            return False


