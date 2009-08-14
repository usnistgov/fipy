class Parallel:
    try:

        try:
            import scipy
        except:
            pass

        from PyTrilinos import Epetra

        procID = Epetra.PyComm().MyPID()
        Nproc = Epetra.PyComm().NumProc()
    except ImportError:
        procID = 0
        Nproc = 1

class Serial:
    procID = 0
    Nproc = 1

parallel = Parallel()
serial = Serial()
import dump
import numerix
import vector
from dimensions.physicalField import PhysicalField
from numerix import *
