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
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this software is not subject to copyright
 # protection and is in the public domain.  FiPy is an experimental
 # system.  NIST assumes no responsibility whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 # 
 # This software can be redistributed and/or modified freely
 # provided that any derivative works bear some notice that they are
 # derived from it, and any modified versions bear some notice that
 # they have been modified.
 # ========================================================================
 #  See the file "license.terms" for information on usage and  redistribution
 #  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #  
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

import cPickle
import os
import gzip

# TODO: add test to show that round trip pickle of mesh doesn't work properly
# FIXME: pickle fails to work properly on numpy 1.1 (run gapFillMesh.py)
def write(data, filename = None, extension = ''):
    """
    Pickle an object and write it to a file. Wrapper for
    `cPickle.dump()`.

    :Parameters:
      - `data`: The object to be pickled.
      - `filename`: The name of the file to place the pickled object. If `filename` is `None`
        then a temporary file will be used and the file object and file name will be returned as a tuple
      - `extension`: Used if filename is not given.

    Test to check pickling and unpickling.

        >>> from fipy.meshes.grid1D import Grid1D
        >>> f, tempfile = write(Grid1D(nx = 2))
        >>> mesh = read(tempfile, f)
        >>> print mesh.getNumberOfCells()
        2
        
    """
    if filename is None:
        import tempfile
        (f, _filename) =  tempfile.mkstemp(extension)
    else:
        _filename = filename
    fileStream = gzip.GzipFile(filename = _filename, mode = 'w', fileobj = None)
    cPickle.dump(data, fileStream, 0)
    fileStream.close()
    if filename is None:
        return (f, _filename)

def read(filename, fileobject = None):
    """
    Read a pickled object from a file. Returns the unpickled object.
    Wrapper for `cPickle.load()`.

    :Parameters:
      - `filename`: The name of the file to unpickle the object from.
      - `fileobject`: Used to remove temporary files
      
    """
    fileStream = gzip.GzipFile(filename = filename, mode = 'r', fileobj = None)
    data = cPickle.load(fileStream)
    fileStream.close()
    if fileobject is not None:
        os.close(fileobject)
        os.remove(filename)
    return data

def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__": 
    _test()     
