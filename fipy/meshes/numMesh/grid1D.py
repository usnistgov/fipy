#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "grid1D.py"
 #                                    created: 11/10/03 {3:30:42 PM} 
 #                                last update: 4/4/05 {11:57:43 AM} 
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

"""1D rectangular Mesh
"""

import Numeric

from fipy.meshes.numMesh.mesh1D import Mesh1D
from fipy.meshes.numMesh.face import Face
from fipy.tools import vector
from fipy.tools.dimensions.physicalField import PhysicalField

class Grid1D(Mesh1D):
    """
    Creates a 1D grid mesh. Vertices and cells are numbered 
    in the usual way.
    
        >>> mesh = Grid1D(nx = 3)
        >>> print mesh.getCellCenters()
        [[ 0.5,]
         [ 1.5,]
         [ 2.5,]] 1
         
        >>> mesh = Grid1D(dx = (1, 2, 3))
        >>> print mesh.getCellCenters()
        [[ 0.5,]
         [ 2. ,]
         [ 4.5,]] 1
         
        >>> mesh = Grid1D(nx = 2, dx = (1, 2, 3))
        Traceback (most recent call last):
        ...
        IndexError: nx != len(dx)

    """
    def __init__(self, dx = 1., nx = None):
        self.dx = PhysicalField(value = dx)
        scale = PhysicalField(value = 1, unit = self.dx.getUnit())
        self.dx /= scale
        
        self.nx = self._calcNumPts(d = self.dx, n = nx)
        
        self.numberOfVertices = self.nx + 1
        
        vertices = self._createVertices()
        faces = self._createFaces()
        cells = self._createCells()
        Mesh1D.__init__(self, vertices, faces, cells)
        
        self.setScale(value = scale)
        
    def _createVertices(self):
        x = self._calcVertexCoordinates(self.dx, self.nx)
        
        return Numeric.transpose(Numeric.array((x,)))
    
    def _createFaces(self):
        return Numeric.arange(self.numberOfVertices)

    def _createCells(self):
        """
        cells = (f1, f2) going left to right.
        f1 etc. refer to the faces
        """
        self.numberOfFaces = self.nx + 1
        f1 = Numeric.arange(self.nx)
        f2 = f1 + 1
        return Numeric.transpose(Numeric.array((f1, f2)))

    def getDim(self):
        return 1
        
    def getFacesLeft(self):
        """Return face on left boundary of Grid1D as list.
        """
        return [Face(self, 0)]
        
    def getFacesRight(self):
        """Return face on right boundary of Grid1D as list.
        """
        return [Face(self, self.numberOfFaces - 1)]
        
    def getScale(self):
        return self.scale['length']
        
    def getPhysicalShape(self):
        """Return physical dimensions of Grid1D.
        """
        return PhysicalField(value = (self.nx * self.dx * self.getScale(),))

    def _getMeshSpacing(self):
        return Numeric.array((self.dx,))
    
    def getShape(self):
        return (self.nx,)
    
## pickling

    def __getstate__(self):
        dict = {
            'dx' : self.dx,            
            'nx' : self.nx}
        return dict

    def __setstate__(self, dict):
        self.__init__(dx = dict['dx'], nx = dict['nx'])

def _test():
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
