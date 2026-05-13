from fipy.tools import numerix
from fipy.tools.comms.commWrapper import CommWrapper

__all__ = ["DummyComm"]

class DummyComm(CommWrapper):
    @property
    def procID(self):
        return 0
        
    @property
    def Nproc(self):
        return 1
