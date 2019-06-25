from __future__ import unicode_literals
from fipy.tools import numerix
from fipy.tools.comms.serialCommWrapper import SerialCommWrapper

__all__ = ["DummyComm"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class DummyComm(SerialCommWrapper):
    def __init__(self):
        pass

    def Barrier(self):
        pass

    def sum(self, a, axis=None):
        return a.sum(axis=axis)

    def MaxAll(self, vec):
        return max(numerix.array(vec))

    def MinAll(self, vec):
        return min(numerix.array(vec))

    def __setstate__(self, dict):
        self.__init__()
