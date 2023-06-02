from builtins import object
__docformat__ = 'restructuredtext'

from fipy.tools import numerix

__all__ = ["CommWrapper"]

class CommWrapper(object):
    """MPI Communicator wrapper
    
    Encapsulates capabilities needed for possibly parallel operations. 
    Some capabilities are not parallel.
    
    """
    
    def __init__(self):
        pass
        
    def __repr__(self):
        return "%s()" % self.__class__.__name__

    @property
    def procID(self):
        raise NotImplementedError
        
    @property
    def Nproc(self):
        raise NotImplementedError
        
    def Barrier(self):
        pass

    def all(self, a, axis=None):
        return a.all(axis=axis)

    def any(self, a, axis=None):
        return a.any(axis=axis)

    def allclose(self, a, b, rtol=1.e-5, atol=1.e-8):
        return numerix.allclose(a, b, rtol=rtol, atol=atol)

    def allequal(self, a, b):
        return numerix.allequal(a, b)

    def bcast(self, obj, root=0):
        return obj

    def allgather(self, obj):
        return obj

    def sum(self, a, axis=None):
        return a.sum(axis=axis)

    def __getstate__(self):
        return {'dummy': 0}

    def __setstate__(self, dict):
        self.__init__()
        
    def Norm2(self, vec):
        return numerix.L2norm(vec)

    @staticmethod
    def _tolist(vec):
        try:
            return list(vec)
        except TypeError:
            return [vec]

    def MaxAll(self, vec):
        return max(self._tolist(vec))
        
    def MinAll(self, vec):
        return min(self._tolist(vec))
