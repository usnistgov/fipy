#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "tri2D.py"
 #                                    created: 07/07/04 {4:28:00 PM} 
 #                                last update: 10/22/04 {4:21:13 PM} 
 #  Author: Alexander Mont <alexander.mont@nist.gov>
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
 #  2004-07-07 ADM 1.0 original
 # ###################################################################
 ##


"""

2-dimensional triangular mesh.

This class creates a mesh made out of triangles. It does this by starting with a standard Cartesian mesh (Grid2D) and dividing each cell in that mesh (hereafter referred to as a 'box') into four equal partswith the dividing lines being the diagonals. A Tri2D mesh is constructed using the following keyword arguments:

dx, dy - The X and Y dimensions of each box. Note that if dx \ne dy, the line segments connecting the cell centers will not be orthogonal to the faces.

nx, ny - The number of boxes in the X direction and the Y direction. The total number of boxes will be equal to nx * ny, and the total number of cells will be equal to 4 * nx * ny.

The faces, cells, and vertices are numbered as follows:

Faces - Horizontal faces are numbered first, then vertical faces, then diagonal faces on the lower left of the boxes, then diagonal faces on the llower right of boxes, then diagonal faces on the upper left of boxes, then diagonal faces on the upper right of boxes.

Vertices - Vertices on the corners of boxes are numbered first, then vertices on the box centers.

Cells - Cells on the right of boxes are numbered first, then cells on the top of boxes, then cells on the left of boxes, then cells on the bottom of boxes.

Within each of the 'sub-categories' in the above, the vertices, cells and faces are numbered in the usual way.



Test cases:

   >>> testmesh = Tri2D(dx = 0.5, dy = 0.5, nx = 2, ny = 2)
   >>> list = testmesh.createVertices().tolist()
   >>> print list
   [[0.0, 0.0], [0.5, 0.0], [1.0, 0.0], [0.0, 0.5], [0.5, 0.5], [1.0, 0.5], [0.0, 1.0], [0.5, 1.0], [1.0, 1.0], [0.25, 0.25], [0.75, 0.25], [0.25, 0.75], [0.75, 0.75]]

   >>> testmesh = Tri2D(dx = 0.5, dy = 0.5, nx = 2, ny = 2)
   >>> list = testmesh.createFaces().tolist()
   >>> print list
   [[1, 0], [2, 1], [3, 4], [4, 5], [6, 7], [7, 8], [0, 3], [4, 1], [5, 2], [3, 6], [7, 4], [8, 5], [9, 0], [10, 1], [11, 3], [12, 4], [1, 9], [2, 10], [4, 11], [5, 12], [9, 3], [10, 4], [11, 6], [12, 7], [9, 4], [10, 5], [11, 7], [12, 8]]

   >>> testmesh = Tri2D(dx = 0.5, dy = 0.5, nx = 2, ny = 2)
   >>> list = testmesh.createCells().tolist()
   >>> print list
   [[7, 24, 16], [8, 25, 17], [10, 26, 18], [11, 27, 19], [2, 20, 24], [3, 21, 25], [4, 22, 26], [5, 23, 27], [6, 12, 20], [7, 13, 21], [9, 14, 22], [10, 15, 23], [0, 16, 12], [1, 17, 13], [2, 18, 14], [3, 19, 15]]

"""

__docformat__ = "restructuredtext"

import Numeric

from fipy.meshes.numMesh.mesh2D import Mesh2D
from fipy.meshes.numMesh.face import Face
import fipy.tools.vector as vector
from fipy.tools.dimensions.physicalField import PhysicalField

class Tri2D(Mesh2D):
    """
    Creates a 2D triangular mesh with horizontal faces numbered
    first then vertical faces, then diagonal faces. Vertices are numbered starting with the
    vertices at the corners of boxes and then the vertices at the centers of boxes. Vertices and cells are numbered 
    in the usual way.
    """
    def __init__(self, dx = 1., dy = 1., nx = None, ny = 1):
        self.nx = nx
        self.ny = ny
	
	self.dx = PhysicalField(value = dx)
	scale = PhysicalField(value = 1, unit = self.dx.getUnit())
	self.dx /= scale
	
	self.dy = PhysicalField(value = dy)
	if self.dy.getUnit().isDimensionless():
	    self.dy = dy
	else:
	    self.dy /= scale
	
        self.numberOfCornerVertices = (self.nx + 1) * (self. ny + 1)
        self.numberOfCenterVertices = self.nx * self.ny
        self.numberOfTotalVertices = self.numberOfCornerVertices + self.numberOfCenterVertices
	
	vertices = self.createVertices()
        faces = self.createFaces()
        cells = self.createCells()
        cells = Numeric.sort(cells)
        Mesh2D.__init__(self, vertices, faces, cells)
	self.setScale(value = scale)
        
    def createVertices(self):
        
        x = Numeric.arange(self.nx + 1) * self.dx
        y = Numeric.arange(self.ny + 1) * self.dy
        x = Numeric.resize(x, (self.numberOfCornerVertices,))
        y = Numeric.repeat(y, self.nx + 1)
        boxCorners = Numeric.transpose(Numeric.array((x, y)))
        x = Numeric.arange(0.5, self.nx + 0.5) * self.dx
        y = Numeric.arange(0.5, self.ny + 0.5) * self.dy
        x = Numeric.resize(x, (self.numberOfCenterVertices,))
        y = Numeric.repeat(y, self.nx)
        boxCenters = Numeric.transpose(Numeric.array((x, y)))
        return Numeric.concatenate((boxCorners, boxCenters))
    
    def createFaces(self):
        """
        v1, v2 refer to the cells.
        Horizontel faces are first
        """
        v1 = Numeric.arange(self.numberOfCornerVertices)
        v2 = v1 + 1
        horizontalFaces = vector.prune(Numeric.transpose(Numeric.array((v1, v2))), self.nx + 1, self.nx)
        v1 = Numeric.arange(self.numberOfCornerVertices - (self.nx + 1))
        v2 = v1 + self.nx + 1
        verticalFaces =  Numeric.transpose(Numeric.array((v1, v2)))

        ## reverse some of the face orientations to obtain the correct normals

        tmp = horizontalFaces.copy()
        horizontalFaces[:self.nx, 0] = tmp[:self.nx, 1]
        horizontalFaces[:self.nx, 1] = tmp[:self.nx, 0]

        tmp = verticalFaces.copy()
        verticalFaces[:, 0] = tmp[:, 1]
        verticalFaces[:, 1] = tmp[:, 0]
        verticalFaces[::(self.nx + 1), 0] = tmp[::(self.nx + 1), 0]
        verticalFaces[::(self.nx + 1), 1] = tmp[::(self.nx + 1), 1]

        ## do the center ones now
        
        cellCenters = Numeric.arange(self.numberOfCornerVertices, self.numberOfTotalVertices)
        lowerLefts = vector.prune(Numeric.arange(self.numberOfCornerVertices - (self.nx + 1)), self.nx + 1, self.nx)
        lowerRights = lowerLefts + 1
        upperLefts = lowerLefts + self.nx + 1
        upperRights = lowerLefts + self.nx + 2
        lowerLeftFaces = Numeric.transpose(Numeric.array((cellCenters, lowerLefts)))
        lowerRightFaces = Numeric.transpose(Numeric.array((lowerRights, cellCenters)))
        upperLeftFaces = Numeric.transpose(Numeric.array((cellCenters, upperLefts)))
        upperRightFaces = Numeric.transpose(Numeric.array((cellCenters, upperRights)))
        return Numeric.concatenate((horizontalFaces, verticalFaces, lowerLeftFaces, lowerRightFaces, upperLeftFaces, upperRightFaces))

    def createCells(self):
        """
        cells = (f1, f2, f3, f4) going anticlockwise.
        f1 etx refer to the faces
        """
        self.numberOfHorizontalFaces = self.nx * (self.ny + 1)
        self.numberOfVerticalFaces =  self.ny * (self.nx + 1)
        self.numberOfEachDiagonalFaces = self.nx * self.ny
        bottomFaces = Numeric.arange(0, self.numberOfHorizontalFaces - self.nx)
        topFaces = Numeric.arange(self.nx, self.numberOfHorizontalFaces)
        leftFaces = vector.prune(Numeric.arange(self.numberOfHorizontalFaces, self.numberOfHorizontalFaces + self.numberOfVerticalFaces), self.nx + 1, self.nx)
        rightFaces = vector.prune(Numeric.arange(self.numberOfHorizontalFaces, self.numberOfHorizontalFaces + self.numberOfVerticalFaces), self.nx + 1, 0)
        lowerLeftDiagonalFaces = Numeric.arange(self.numberOfHorizontalFaces + self.numberOfVerticalFaces, self.numberOfHorizontalFaces + self.numberOfVerticalFaces + self.numberOfEachDiagonalFaces)
        lowerRightDiagonalFaces = lowerLeftDiagonalFaces + self.numberOfEachDiagonalFaces
        upperLeftDiagonalFaces = lowerRightDiagonalFaces + self.numberOfEachDiagonalFaces
        upperRightDiagonalFaces = upperLeftDiagonalFaces + self.numberOfEachDiagonalFaces
        ##faces in arrays, now get the cells
        bottomOfBoxCells = Numeric.transpose(Numeric.array([bottomFaces, lowerRightDiagonalFaces, lowerLeftDiagonalFaces]))
        rightOfBoxCells = Numeric.transpose(Numeric.array([rightFaces, upperRightDiagonalFaces, lowerRightDiagonalFaces]))
        topOfBoxCells = Numeric.transpose(Numeric.array([topFaces, upperLeftDiagonalFaces, upperRightDiagonalFaces]))
        leftOfBoxCells = Numeric.transpose(Numeric.array([leftFaces, lowerLeftDiagonalFaces, upperLeftDiagonalFaces]))
        return Numeric.concatenate((rightOfBoxCells, topOfBoxCells, leftOfBoxCells, bottomOfBoxCells))

    def getFacesLeft(self):
	"""Return list of faces on left boundary of Grid2D.
	"""
	return [Face(self, id) for id in Numeric.arange(self.numberOfHorizontalFaces, self.numberOfHorizontalFaces + self.numberOfVerticalFaces, self.nx + 1)]
	
    def getFacesRight(self):
	"""Return list of faces on right boundary of Grid2D.
	"""
	return [Face(self, id) for id in  Numeric.arange(self.numberOfHorizontalFaces + self.nx, self.numberOfHorizontalFaces + self.numberOfVerticalFaces, self.nx + 1)]
	
    def getFacesTop(self):
	"""Return list of faces on top boundary of Grid2D.
	"""
	return [Face(self, id) for id in Numeric.arange(self.numberOfHorizontalFaces - self.nx, self.numberOfHorizontalFaces)]
	
    def getFacesBottom(self):
	"""Return list of faces on bottom boundary of Grid2D.
	"""
	return [Face(self, id) for id in Numeric.arange(self.nx)]
        
    def getScale(self):
	return self.scale['length']
	
    def getPhysicalShape(self):
	"""Return physical dimensions of Grid2D.
	"""
	return PhysicalField(value = (self.nx * self.dx * self.getScale(), self.ny * self.dy * self.getScale()))



    def getMeshSpacing(self):
	return Numeric.array((self.dx,self.dy))
    
    def getShape(self):
        return (self.nx, self.ny)
    
## pickling

    def __getstate__(self):
        dict = {
            'dx' : self.dx,            
            'dy' : self.dy,
            'nx' : self.nx,
            'ny' : self.ny}
        return dict

    def __setstate__(self, dict):
        self.__init__(dx = dict['dx'], dy = dict['dy'], nx = dict['nx'], ny = dict['ny'])

## test test test
        
def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__": 
    _test() 







