"""This module provides a generic, high-level interface for
creating shared temporary files.  All of the interfaces
provided by this module can be used without fear of race conditions.
"""

from tempfile import NamedTemporaryFile, template, _TemporaryFileWrapper
from future import utils

from fipy.tools import parallelComm

__all__ = ["SharedTemporaryFile"]

# Silence linting errors that come from following conventions of
# tempfile package.
# pylint: disable=invalid-name
# pylint: disable=redefined-builtin
# pylint: disable=too-many-arguments

if utils.PY3:
    import io as _io

    def SharedTemporaryFile(
        mode="w+b",
        buffering=-1,
        encoding=None,
        newline=None,
        suffix="",
        prefix=template,
        dir=None,
        delete=True,
        communicator=parallelComm,
    ):
        """Create a temporary file shared by all MPI ranks.

        The file is created as `NamedTemporaryFile` would do it.
        The name of the returned file-like object is accessible as its
        ``name`` attribute.  The file will be automatically deleted when it
        is closed unless the `delete` argument is set to False.

        >>> from fipy.tools import SharedTemporaryFile, parallelComm
        >>> with SharedTemporaryFile(mode='w+', suffix=".tmp") as tmpFile:
        ...     # write on processor 0
        ...     if parallelComm.procID == 0:
        ...         _ = tmpFile.write("shared text")
        ...
        ...     parallelComm.Barrier()
        ...
        ...     # read on all processors
        ...     _ = tmpFile.seek(0)
        ...     txt = tmpFile.read()
        >>> print(txt)
        shared text

        Parameters
        ----------
        prefix, suffix, dir : str
            As for mkstemp
        mode : str
            The mode argument to io.open (default "w+b")
        buffering : int
            The buffer size argument to io.open (default -1)
        encoding : str or None
            The encoding argument to io.open (default None)
        newline : str or None
            The newline argument to io.open (default None)
        delete : bool
            Whether the file is deleted on close (default True)
        communicator : ~fipy.tools.comms.commWrapper.CommWrapper
            MPI communicator describing ranks to share with.  A duck-typed
            object with `procID` and `Nproc` attributes is sufficient.

        Returns
        -------
        file-like object

        See Also
        --------
        tempfile.NamedTemporaryFile, tempfile.mkstemp, io.open

        """

        if communicator.procID == 0:
            file = NamedTemporaryFile(
                mode=mode,
                buffering=buffering,
                encoding=encoding,
                newline=newline,
                suffix=suffix,
                prefix=prefix,
                dir=dir,
                delete=delete,
            )
            fname = file.name
        else:
            fname = None

        fname = communicator.bcast(fname)

        if communicator.procID != 0:
            file = _io.open(
                fname, mode, buffering=buffering, newline=newline, encoding=encoding
            )
            # let procID 0 handle delete
            file = _TemporaryFileWrapper(file, fname, delete=False)

        return file


else:

    def SharedTemporaryFile(
        mode="w+b",
        bufsize=-1,
        suffix="",
        prefix=template,
        dir=None,
        delete=True,
        communicator=parallelComm,
    ):
        """Create a temporary file shared by all MPI ranks.

        The file is created as `NamedTemporaryFile` would do it.
        The name of the returned file-like object is accessible as its
        ``name`` attribute.  The file will be automatically deleted when it
        is closed unless the `delete` argument is set to False.

        >>> from fipy.tools import SharedTemporaryFile, parallelComm
        >>> with SharedTemporaryFile(mode='w+', suffix=".tmp") as tmpFile:
        ...     # write on processor 0
        ...     if parallelComm.procID == 0:
        ...         _ = tmpFile.write("shared text")
        ...
        ...     parallelComm.Barrier()
        ...
        ...     # read on all processors
        ...     _ = tmpFile.seek(0)
        ...     txt = tmpFile.read()
        >>> print(txt)
        shared text

        Parameters
        ----------
        prefix, suffix, dir : str
            As for mkstemp
        mode : str
            The mode argument to os.fdopen (default "w+b")
        bufsize : int
            The buffer size argument to os.fdopen (default -1)
        delete : bool
            Whether the file is deleted on close (default True)
        communicator : ~fipy.tools.comms.commWrapper.CommWrapper
            MPI communicator describing ranks to share with.  A duck-typed
            object with `procID` and `Nproc` attributes is sufficient

        Returns
        -------
        file-like object

        See Also
        --------
        tempfile.NamedTemporaryFile, tempfile.mkstemp, os.fdopen

        """

        if communicator.procID == 0:
            file = NamedTemporaryFile(
                mode=mode,
                bufsize=bufsize,
                suffix=suffix,
                prefix=prefix,
                dir=dir,
                delete=delete,
            )
            fname = file.name
        else:
            fname = None

        fname = communicator.bcast(fname)

        if communicator.procID != 0:
            file = open(fname, mode, bufsize)
            # let procID 0 handle delete
            file = _TemporaryFileWrapper(file, fname, delete=False)

        return file


def _test():
    import fipy.tests.doctestPlus

    return fipy.tests.doctestPlus.testmod()


if __name__ == "__main__":
    _test()
