#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "dump.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # of Standards and Technology, an agency of the Federal Government.
 # Pursuant to title 17 section 105 of the United States Code,
 # United States Code this software is not subject to copyright
 # protection, and this software is considered to be in the public domain.
 # FiPy is an experimental system.
 # NIST assumes no responsibility whatsoever for its use by whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 #
 # To the extent that NIST may hold copyright in countries other than the
 # United States, you are hereby granted the non-exclusive irrevocable and
 # unconditional right to print, publish, prepare derivative works and
 # distribute this software, in any medium, or authorize others to do so on
 # your behalf, on a royalty-free basis throughout the world.
 #
 # You may improve, modify, and create derivative works of the software or
 # any portion of the software, and you may copy and distribute such
 # modifications or works.  Modified works should carry a notice stating
 # that you changed the software and should note the date and nature of any
 # such change.  Please explicitly acknowledge the National Institute of
 # Standards and Technology as the original source.
 #
 # This software can be redistributed and/or modified freely provided that
 # any derivative works bear some notice that they are derived from it, and
 # any modified versions bear some notice that they have been modified.
 # ========================================================================
 #
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

import cPickle
import os
import sys
import gzip

from fipy.tools import parallelComm

__all__ = ["write", "read"]

# TODO: add test to show that round trip pickle of mesh doesn't work properly
# FIXME: pickle fails to work properly on numpy 1.1 (run gapFillMesh.py)
def write(data, filename = None, extension = '', communicator=parallelComm):
    """
    Pickle an object and write it to a file. Wrapper for
    `cPickle.dump()`.

    :Parameters:
      - `data`: The object to be pickled.
      - `filename`: The name of the file to place the pickled object. If `filename` is `None`
        then a temporary file will be used and the file object and file name will be returned as a tuple
      - `extension`: Used if filename is not given.
      - `communicator`: Object with `procID` and `Nproc` attributes.

    Test to check pickling and unpickling.

        >>> from fipy.meshes import Grid1D
        >>> old = Grid1D(nx = 2)
        >>> f, tempfile = write(old)
        >>> new = read(tempfile, f)
        >>> print old.numberOfCells == new.numberOfCells
        True

    """
    if communicator.procID == 0:
        if filename is None:
            import tempfile
            (f, _filename) =  tempfile.mkstemp(extension)
        else:
            (f, _filename) = (None, filename)
        fileStream = gzip.GzipFile(filename = _filename, mode = 'w', fileobj = None)
    else:
        fileStream = open(os.devnull, mode='w')
        (f, _filename) = (None, os.devnull)

    cPickle.dump(data, fileStream, 0)
    fileStream.close()

    if filename is None:
        return (f, _filename)

def read(filename, fileobject=None, communicator=parallelComm, mesh_unmangle=False):
    """
    Read a pickled object from a file. Returns the unpickled object.
    Wrapper for `cPickle.load()`.

    :Parameters:
      - `filename`: The name of the file to unpickle the object from.
      - `fileobject`: Used to remove temporary files
      - `communicator`: Object with `procID` and `Nproc` attributes.
      - `mesh_unmangle`: Correct improper pickling of non-uniform meshes (ticket:243)

    """
    if communicator.procID == 0:
        fileStream = gzip.GzipFile(filename = filename, mode = 'r', fileobj = None)
        data = fileStream.read()
        fileStream.close()
        if fileobject is not None:
            os.close(fileobject)
            os.remove(filename)
    else:
        data = None

    if communicator.Nproc > 1:
        data = communicator.bcast(data, root=0)

    if sys.version_info < (3,0):
        import StringIO
        f = StringIO.StringIO(data)
    else:
        import io
        f = io.BytesIO(data)

    unpickler = cPickle.Unpickler(f)

    if mesh_unmangle:
        def find_class(module, name):
            __import__(module)
            mod = sys.modules[module]
            klass = getattr(mod, name)

            from fipy import meshes
            import types

            if isinstance(klass, types.ClassType) and issubclass(klass, meshes.mesh.Mesh):
                class UnmangledMesh(klass):
                    def __setstate__(self, dict):
                        if ('cellFaceIDs' in dict
                            and 'faceVertexIDs' in dict):

                            dict = dict.copy()
                            for key in ('cellFaceIDs', 'faceVertexIDs'):
                                arr = dict[key]
                                arr.data[:] = arr.transpose().flatten().reshape(arr.shape)

                        klass.__setstate__(self, dict)

                return UnmangledMesh
            else:
                return klass

        unpickler.find_global = find_class

    return unpickler.load()

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
