


from fipy.tools.comms.commWrapper import CommWrapper

__all__ = ["SerialCommWrapper"]

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
