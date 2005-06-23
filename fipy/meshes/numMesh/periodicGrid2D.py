#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "periodicGrid2D.py"
 #                                    created: 11/10/03 {3:30:42 PM} 
 #                                last update: 4/21/05 {4:49:30 PM} 
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
 #  
 #  Description: 
 # 
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2003-11-10 JEG 1.0 original
 # ###################################################################
 ##

"""
2D periodic rectangular Mesh
"""
__docformat__ = 'restructuredtext'

from fipy.meshes.numMesh.grid2D import Grid2D

class PeriodicGrid2D(Grid2D):
    """
    Creates a periodic2D grid mesh with horizontal faces numbered
    first and then vertical faces. Vertices and cells are numbered 
    in the usual way.

        >>> mesh = PeriodicGrid2D(dx = 1., dy = 0.5, nx = 2, ny = 2)

        >>> print mesh.getExteriorFaceIDs()
        [ 4, 5, 8,11,]

        >>> print mesh.getFaceCellIDs()
        [[2 ,0 ,]
         [3 ,1 ,]
         [0 ,2 ,]
         [1 ,3 ,]
         [2 ,-- ,]
         [3 ,-- ,]
         [1 ,0 ,]
         [0 ,1 ,]
         [1 ,-- ,]
         [3 ,2 ,]
         [2 ,3 ,]
         [3 ,-- ,]]

        >>> print mesh._getCellDistances()
        [ 0.5 , 0.5 , 0.5 , 0.5 , 0.25, 0.25, 1.  , 1.  , 0.5 , 1.  , 1.  , 0.5 ,] 1
 
        >>> print mesh._getCellFaceIDs()
        [[ 0, 7, 2, 6,]
         [ 1, 6, 3, 7,]
         [ 2,10, 0, 9,]
         [ 3, 9, 1,10,]]

        >>> print mesh._getCellToCellDistances()
        [[ 0.5, 1. , 0.5, 1. ,]
         [ 0.5, 1. , 0.5, 1. ,]
         [ 0.5, 1. , 0.5, 1. ,]
         [ 0.5, 1. , 0.5, 1. ,]] 1

        >>> print mesh._getFaceNormals()
        [[-0., 1.,]
         [-0., 1.,]
         [-0., 1.,]
         [-0., 1.,]
         [-0., 1.,]
         [-0., 1.,]
         [ 1., 0.,]
         [ 1., 0.,]
         [ 1., 0.,]
         [ 1., 0.,]
         [ 1., 0.,]
         [ 1., 0.,]]


    """
    def __init__(self, dx = 1., dy = 1., nx = None, ny = None):
        Grid2D.__init__(self, dx = dx, dy = dy, nx = nx, ny = ny)
        self._connectFaces(self.getFacesLeft(), self.getFacesRight())
        self._connectFaces(self.getFacesBottom(), self.getFacesTop())

def _test():
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
