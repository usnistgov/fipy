from __future__ import unicode_literals

import logging

_log = logging.getLogger(__name__)

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
    # When run in MPI mode, the first Trilinos import makes the "current
    # directory" be the directory with the executable file that's being
    # run.  As best I can tell, this happens in MPI_Init, deep in Trilinos.
    # Possibly because "current directory" not well-defined in MPI between
    # processors?

    # This fix relies on this being the FIRST place to import any Trilinos
    # module.  The only way to import Trilinos things should be to do
    # `from fipy.solvers import *` and have it automatically import Trilinos
    # via this file.

    from PyTrilinos import Epetra

    try:
        import platform
        if platform.dist()[0] == 'debian':
            import PyTrilinos
            if '10.0.4' in PyTrilinos.version():
                # The package mpi4py is a required package if you are using
                # Trilinos on a Debian platform with Trilinos version 10.0.4 due to
                # a Trilinos bug (see <https://github.com/usnistgov/fipy/issues/301>).

                from mpi4py import MPI
    except:
        pass

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
