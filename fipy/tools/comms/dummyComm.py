from __future__ import unicode_literals
from fipy.tools import numerix
from fipy.tools.comms.abstractCommWrapper import AbstractCommWrapper

__all__ = ["DummyComm"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class DummyComm(AbstractCommWrapper):
    @property
    def procID(self):
        return 0
        
    @property
    def Nproc(self):
        return 1
