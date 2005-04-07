#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "dump.py"
 #                                    created: 1/10/04 {10:23:17 AM} 
 #                                last update: 9/3/04 {10:35:28 PM} 
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

def write(data, fileName):
    """
    Pickle an object and write it to a file. Wrapper for
    `cPickle.dump`.

    :Parameters:
      - `data` : The object to be pickled.
      - `fileName` : The name of the file to place the pickled object.

    Test to check pickling and unpickling.

        >>> from fipy.meshes.grid1D import Grid1D
        >>> import tempfile
        >>> f, tempFile = tempfile.mkstemp('.gz')
        >>> write(Grid1D(nx = 2), tempFile)
        >>> mesh = read(tempFile)
        >>> print mesh.getNumberOfCells()
        2
        
    """
    fileStream = gzip.GzipFile(filename = fileName, mode = 'w', fileobj = None)
    cPickle.dump(data, fileStream, 0)
    fileStream.close()

def read(fileName = None):
    """
    Read a pickled object from a file. Returns the unpickled object.
    Wrapper for `cPickle.load`.

    :Parameters:
      - `fileName` : The name of the file to unpickle teh obkect from.

    """
    fileStream = gzip.GzipFile(filename = fileName, mode = 'r', fileobj = None)
    data = cPickle.load(fileStream)
    fileStream.close()
    return data

def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__": 
    _test()     
