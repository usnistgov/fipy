from future import utils
from tempfile import NamedTemporaryFile, template, _TemporaryFileWrapper

from fipy.tools import parallelComm

__all__ = ["SharedTemporaryFile"]

if utils.PY3:
    import io as _io
    
    def SharedTemporaryFile(mode='w+b', buffering=-1, encoding=None,
                            newline=None, suffix="", prefix=template,
                            dir=None, delete=True):
        """Create and return a temporary file shared by all MPI ranks.
        Arguments:
        'prefix', 'suffix', 'dir' -- as for mkstemp.
        'mode' -- the mode argument to io.open (default "w+b").
        'buffering' -- the buffer size argument to io.open (default -1).
        'encoding' -- the encoding argument to io.open (default None)
        'newline' -- the newline argument to io.open (default None)
        'delete' -- whether the file is deleted on close (default True).
        The file is created as mkstemp() would do it.
        Returns an object with a file-like interface; the name of the file
        is accessible as file.name.  The file will be automatically deleted
        when it is closed unless the 'delete' argument is set to False.
        """

        if parallelComm.procID == 0:
            f = NamedTemporaryFile(mode=mode, buffering=buffering, encoding=encoding,
                                   newline=newline, suffix=suffix, prefix=prefix,
                                   dir=dir, delete=delete)
            fname = f.name
        else:
            fname = None

        fname = parallelComm.bcast(fname)
        
        if parallelComm.procID != 0:
            f = _io.open(fname, mode, buffering=buffering,
                         newline=newline, encoding=encoding)
            # let procID 0 handle delete
            f = _TemporaryFileWrapper(f, fname, delete=False)
        
        return f
else:
    def SharedTemporaryFile(mode='w+b', bufsize=-1, suffix="",
                            prefix=template, dir=None, delete=True):
        """Create and return a temporary file shared by all MPI ranks.
        Arguments:
        'prefix', 'suffix', 'dir' -- as for mkstemp.
        'mode' -- the mode argument to os.fdopen (default "w+b").
        'bufsize' -- the buffer size argument to os.fdopen (default -1).
        'delete' -- whether the file is deleted on close (default True).
        The file is created as mkstemp() would do it.
        Returns an object with a file-like interface; the name of the file
        is accessible as its 'name' attribute.  The file will be automatically
        deleted when it is closed unless the 'delete' argument is set to False.
        """

        if parallelComm.procID == 0:
            f = NamedTemporaryFile(mode=mode, bufsize=bufsize, suffix=suffix, 
                                   prefix=prefix, dir=dir, delete=delete)
            fname = f.name
        else:
            fname = None

        fname = parallelComm.bcast(fname)
        
        if parallelComm.procID != 0:
            f = open(fname, mode, bufsize)
            # let procID 0 handle delete
            f = _TemporaryFileWrapper(f, fname, delete=False)
        
        return f           
