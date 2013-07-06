__docformat__ = 'restructuredtext'

from fipy.tools import numerix

__all__ = ["AbstractCommWrapper"]

class AbstractCommWrapper(object):
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
        raise NotImplementedError

    def __getstate__(self):
        return {'dummy': 0}

    def __setstate__(self, dict):
        self.__init__()
        
    def Norm2(self, vec):
        raise NotImplementedError

    def MaxAll(self, vec):
        raise NotImplementedError
        
    def MinAll(self, vec):
        raise NotImplementedError
