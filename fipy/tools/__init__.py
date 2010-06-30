class Parallel(object):
    def __init__(self):
        try:
            try:
                import scipy
            except:
                pass

            from PyTrilinos import Epetra

            self.procID = Epetra.PyComm().MyPID()
            self.Nproc = Epetra.PyComm().NumProc()
        except ImportError:
            self.procID = 0
            self.Nproc = 1
    
    def sumAll(self, pyObj):
        from PyTrilinos import Epetra
        return Epetra.PyComm().SumAll(pyObj)

class Serial(object):
    procID = 0
    Nproc = 1


parallel = Parallel()
serial = Serial()
import dump
import numerix
import vector
from dimensions.physicalField import PhysicalField
from numerix import *
from vitals import Vitals

