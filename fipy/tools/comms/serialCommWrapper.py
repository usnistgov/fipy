from __future__ import unicode_literals
from fipy.tools.comms.commWrapper import CommWrapper

__all__ = ["SerialCommWrapper"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class SerialCommWrapper(CommWrapper):
    @property
    def procID(self):
        return 0

    @property
    def Nproc(self):
        return 1

    def Norm2(self, vec):
        from fipy.tools import numerix
        return numerix.L2norm(numerix.asarray(vec))
