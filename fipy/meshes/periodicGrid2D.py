#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "periodicGrid2D.py"
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
 #  
 # ###################################################################
 ##

"""
2D periodic rectangular Mesh
"""
__docformat__ = 'restructuredtext'

from fipy.meshes.grid2D import Grid2D

class PeriodicGrid2D(Grid2D):
    """
    Creates a periodic2D grid mesh with horizontal faces numbered
    first and then vertical faces. Vertices and cells are numbered 
    in the usual way.

        >>> from fipy import numerix
        >>> from fipy.tools import parallel

        >>> mesh = PeriodicGrid2D(dx = 1., dy = 0.5, nx = 2, ny = 2)
        
        >>> print (parallel.procID > 0 or 
        ...        numerix.allclose(numerix.nonzero(mesh.exteriorFaces)[0],
        ...                         [ 4,  5,  8, 11]))
        True

        >>> print (parallel.procID > 0 or 
        ...        numerix.allclose(mesh.faceCellIDs.filled(-1),
        ...                         [[2, 3, 0, 1, 2, 3, 1, 0, 1, 3, 2, 3],
        ...                          [0, 1, 2, 3, -1, -1, 0, 1, -1, 2, 3, -1]]))
        True

        >>> print (parallel.procID > 0 or 
        ...        numerix.allclose(mesh._cellDistances,
        ...                         [ 0.5, 0.5, 0.5, 0.5, 0.25, 0.25, 1., 1., 0.5, 1., 1., 0.5]))
        True
 
        >>> print (parallel.procID > 0 
        ...        or numerix.allclose(mesh.cellFaceIDs,
        ...                            [[0, 1, 2, 3],
        ...                             [7, 6, 10, 9],
        ...                             [2, 3, 0, 1],
        ...                             [6, 7, 9, 10]]))
        True

        >>> print (parallel.procID > 0 or 
        ...        numerix.allclose(mesh._cellToCellDistances,
        ...                         [[ 0.5, 0.5, 0.5, 0.5],
        ...                          [ 1., 1., 1., 1. ],
        ...                          [ 0.5, 0.5, 0.5, 0.5],
        ...                          [ 1., 1., 1., 1. ]]))
        True

        >>> normals = [[0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1],
        ...            [1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0]]

        >>> print (parallel.procID > 0 
        ...        or numerix.allclose(mesh._faceNormals, normals))
        True

        >>> print (parallel.procID > 0 
        ...        or numerix.allclose(mesh._cellVertexIDs,
        ...                            [[4, 5, 7, 8],
        ...                             [3, 4, 6, 7],
        ...                             [1, 2, 4, 5],
        ...                             [0, 1, 3, 4]]))
        True
    """
    def __init__(self, dx = 1., dy = 1., nx = None, ny = None):
        Grid2D.__init__(self, dx = dx, dy = dy, nx = nx, ny = ny)
        self.nonPeriodicCellVertexIDs \
                = super(PeriodicGrid2D, self)._cellVertexIDs
        self.nonPeriodicOrderedCellVertexIDs \
                = super(PeriodicGrid2D, self)._orderedCellVertexIDs
        from fipy.tools import numerix
        self._connectFaces(numerix.nonzero(self.getFacesLeft()), 
                           numerix.nonzero(self.getFacesRight()))
        self._connectFaces(numerix.nonzero(self.getFacesBottom()), 
                           numerix.nonzero(self.getFacesTop()))

    @property
    def _cellVertexIDs(self):
        return self.nonPeriodicCellVertexIDs

    @property
    def _orderedCellVertexIDs(self):
        return self.nonPeriodicOrderedCellVertexIDs
               
class PeriodicGrid2DLeftRight(PeriodicGrid2D):
    def __init__(self, dx = 1., dy = 1., nx = None, ny = None):
        Grid2D.__init__(self, dx = dx, dy = dy, nx = nx, ny = ny)
        self.nonPeriodicCellVertexIDs = Grid2D._getCellVertexIDs(self)
        self.nonPeriodicOrderedCellVertexIDs = Grid2D._getOrderedCellVertexIDs(self)
        from fipy.tools import numerix
        self._connectFaces(numerix.nonzero(self.getFacesLeft()), 
                           numerix.nonzero(self.getFacesRight()))

class PeriodicGrid2DTopBottom(PeriodicGrid2D):
    def __init__(self, dx = 1., dy = 1., nx = None, ny = None):
        Grid2D.__init__(self, dx = dx, dy = dy, nx = nx, ny = ny)
        self.nonPeriodicCellVertexIDs = Grid2D._getCellVertexIDs(self)
        self.nonPeriodicOrderedCellVertexIDs = Grid2D._getOrderedCellVertexIDs(self)
        from fipy.tools import numerix
        self._connectFaces(numerix.nonzero(self.getFacesBottom()), 
                           numerix.nonzero(self.getFacesTop()))
    
def _test():
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
