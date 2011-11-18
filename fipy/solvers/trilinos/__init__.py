def _dealWithTrilinosImportPathologies():
    ## The import scipy statement is added to allow PyTrilinos to run
    ## without throwing a segmentation fault. This is caused by weird
    ## behavior in scipy and PyTrilinos depending on the order in which
    ## modules are imported

    try:
        import scipy
    except:
        pass

    # The fact that I have to do the following manipulation with the current
    # directory is really, really bad. 
    import os
    current_working_directory_path = os.getcwd()
    from PyTrilinos import ML # Gets around strange Trilinos import-order bugs. 
    os.chdir(current_working_directory_path)
    # When run in MPI mode, the first Trilinos import makes the "current directory"
    # be the directory with the executable file that's being run.  As best I can
    # tell, this happens in MPI_Init, deep in Trilinos. Possibly because "current
    # directory" not well-defined in MPI between processors?

    # This fix relies on this being the FIRST place to import any Trilinos module.
    # The only way to import Trilinos things should be to do "from fipy.solvers
    # import *" and have it automatically import Trilinos via this file.

_dealWithTrilinosImportPathologies()

from fipy.solvers.trilinos.preconditioners import *

from fipy.solvers.trilinos.linearCGSSolver import *
from fipy.solvers.trilinos.linearPCGSolver import *
from fipy.solvers.trilinos.linearGMRESSolver import *
from fipy.solvers.trilinos.linearLUSolver import *
from fipy.solvers.trilinos.linearBicgstabSolver import *

from fipy.solvers.trilinos.trilinosMLTest import *

DefaultSolver = LinearGMRESSolver
DefaultAsymmetricSolver = LinearGMRESSolver
DummySolver = DefaultSolver
GeneralSolver = DefaultSolver

__all__ = ["DefaultSolver",
           "DummySolver",
           "DefaultAsymmetricSolver",
           "GeneralSolver"]
           
__all__.extend(linearCGSSolver.__all__)
__all__.extend(linearPCGSolver.__all__)
__all__.extend(linearGMRESSolver.__all__)
__all__.extend(linearLUSolver.__all__)
__all__.extend(linearBicgstabSolver.__all__)
__all__.extend(trilinosMLTest.__all__)
__all__.extend(preconditioners.__all__)
