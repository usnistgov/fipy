#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "tri2D.py"
 #                                    created: 07/07/04 {4:28:00 PM} 
 #                                last update: 3/5/06 {8:15:53 AM} 
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

__docformat__ = "restructuredtext"

import Numeric

from fipy.meshes.numMesh.mesh2D import Mesh2D
from fipy.meshes.meshIterator import FaceIterator
from fipy.tools import vector
from fipy.tools.dimensions.physicalField import PhysicalField

class Tri2D(Mesh2D):
    """
    This class creates a mesh made out of triangles.  It does this by
    starting with a standard Cartesian mesh (`Grid2D`) and dividing each cell
    in that mesh (hereafter referred to as a 'box') into four equal
    parts with the dividing lines being the diagonals.
    """
    
    def __init__(self, dx = 1., dy = 1., nx = 1, ny = 1):
        """
        Creates a 2D triangular mesh with horizontal faces numbered first then
        vertical faces, then diagonal faces.  Vertices are numbered starting
        with the vertices at the corners of boxes and then the vertices at the
        centers of boxes.  Cells on the right of boxes are numbered first, then
        cells on the top of boxes, then cells on the left of boxes, then cells
        on the bottom of boxes.  Within each of the 'sub-categories' in the
        above, the vertices, cells and faces are numbered in the usual way.
        
        :Parameters:
          - `dx, dy`: The X and Y dimensions of each 'box'. 
            If `dx` <> `dy`, the line segments connecting the cell 
            centers will not be orthogonal to the faces.
          - `nx, ny`: The number of boxes in the X direction and the Y direction. 
            The total number of boxes will be equal to `nx * ny`, and the total 
            number of cells will be equal to `4 * nx * ny`.
        """
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
	
	vertices = self._createVertices()
        faces = self._createFaces()
        cells = self._createCells()
        cells = Numeric.sort(cells)
        Mesh2D.__init__(self, vertices, faces, cells)
	self.setScale(value = scale)
        
    def _createVertices(self):
        
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
    
    def _createFaces(self):
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

    def _createCells(self):
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
	return FaceIterator(mesh=self,
                            ids=Numeric.arange(self.numberOfHorizontalFaces, 
                                               self.numberOfHorizontalFaces + self.numberOfVerticalFaces, 
                                               self.nx + 1))
	
    def getFacesRight(self):
	"""Return list of faces on right boundary of Grid2D.
	"""
	return FaceIterator(mesh=self,
                            ids=Numeric.arange(self.numberOfHorizontalFaces + self.nx, 
                                               self.numberOfHorizontalFaces + self.numberOfVerticalFaces, 
                                               self.nx + 1))
	
    def getFacesTop(self):
	"""Return list of faces on top boundary of Grid2D.
	"""
	return FaceIterator(mesh=self, 
                            ids=Numeric.arange(self.numberOfHorizontalFaces - self.nx, 
                                               self.numberOfHorizontalFaces))
	
    def getFacesBottom(self):
	"""Return list of faces on bottom boundary of Grid2D.
	"""
	return FaceIterator(mesh=self, 
                            ids=Numeric.arange(self.nx))
        
    def getScale(self):
	return self.scale['length']
	
    def getPhysicalShape(self):
	"""Return physical dimensions of Grid2D.
	"""
	return PhysicalField(value = (self.nx * self.dx * self.getScale(), self.ny * self.dy * self.getScale()))



    def _getMeshSpacing(self):
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

        
    def _test(self):
        """
        These tests are not useful as documentation, but are here to ensure
        everything works as expected.

            >>> dx = 0.5
            >>> dy = 2.
            >>> nx = 3
            >>> ny = 2
            
            >>> mesh = Tri2D(nx = nx, ny = ny, dx = dx, dy = dy)     
            
            >>> vertices = Numeric.array(((0.0, 0.0), (0.5, 0.0), (1.0, 0.0), (1.5, 0.0),
            ...                           (0.0, 2.0), (0.5, 2.0), (1.0, 2.0), (1.5, 2.0),
            ...                           (0.0, 4.0), (0.5, 4.0), (1.0, 4.0), (1.5, 4.0),
            ...                           (0.25, 1.0), (0.75, 1.0), (1.25, 1.0),
            ...                           (0.25, 3.0), (0.75, 3.0), (1.25, 3.0)))
            
            >>> from fipy.tools import numerix
            >>> numerix.allequal(vertices, mesh._createVertices())
            1
        
            >>> faces = Numeric.array(((1, 0), (2, 1), (3, 2),
            ...                        (4, 5), (5, 6), (6, 7),
            ...                        (8, 9), (9, 10), (10, 11),
            ...                        (0, 4), (5, 1), (6, 2), (7, 3),
            ...                        (4, 8), (9, 5), (10, 6), (11, 7),
            ...                        (12, 0), (13, 1), (14, 2), (15, 4), (16, 5), (17, 6),
            ...                        (1, 12), (2, 13), (3, 14), (5, 15), (6, 16), (7, 17),
            ...                        (12, 4), (13, 5), (14, 6), (15, 8), (16, 9), (17, 10),
            ...                        (12, 5), (13, 6), (14, 7), (15, 9), (16, 10), (17, 11)))
            >>> numerix.allequal(faces, mesh._createFaces())
            1

            >>> cells = Numeric.array(((10, 35, 23), (11, 36, 24), (12, 37, 25), (14, 38, 26), (15, 39, 27), (16, 40, 28),
            ...                        (3, 29, 35), (4, 30, 36), (5, 31, 37), (6, 32, 38), (7, 33, 39), (8, 34, 40),
            ...                        (9, 17, 29), (10, 18, 30), (11, 19, 31), (13, 20, 32), (14, 21, 33), (15, 22, 34),
            ...                        (0, 23, 17), (1, 24, 18), (2, 25, 19), (3, 26, 20), (4, 27, 21), (5, 28, 22)))
            >>> numerix.allequal(cells, mesh._createCells())
            1

            >>> externalFaces = Numeric.array((0, 1, 2, 6, 7, 8, 9 , 12, 13, 16))
            >>> numerix.allequal(externalFaces, mesh.getExteriorFaces())
            1

            >>> internalFaces = Numeric.array((3, 4, 5, 10, 11, 14, 15, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40))
            >>> numerix.allequal(internalFaces, mesh.getInteriorFaces())
            1

            >>> import MA
            >>> faceCellIds = MA.masked_values(((18, -1), (19, -1), (20, -1), (6, 21), (7, 22), (8, 23), (9, -1), (10, -1), (11, -1),
            ...                                 (12, -1), (0, 13), (1, 14), (2, -1), (15, -1), (3, 16), (4, 17), (5, -1),
            ...                                 (12, 18), (13, 19), (14, 20), (15, 21), (16, 22), (17, 23),
            ...                                 (0, 18), (1, 19), (2, 20), (3, 21), (4, 22), (5, 23),
            ...                                 (6, 12), (7, 13), (8, 14), (9, 15), (10, 16), (11, 17),
            ...                                 (0, 6), (1, 7), (2, 8), (3, 9), (4, 10), (5, 11)), -1)
            >>> numerix.allequal(faceCellIds, mesh.getFaceCellIDs())
            1
            
            >>> d = (Numeric.sqrt((dx*dx)+(dy*dy))) / 2.0 ## length of diagonal edges  
            >>> faceAreas = Numeric.array((0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
            ...                            2, 2, 2, 2, 2, 2, 2, 2,
            ...                            d, d, d, d, d, d, d, d, d, d, d, d, d, d, d, d, d, d, d, d, d, d, d, d))
            >>> numerix.allclose(faceAreas, mesh._getFaceAreas(), atol = 1e-10, rtol = 1e-10)
            1
            
            >>> faceCoords = Numeric.take(vertices, faces)
            >>> faceCenters = (faceCoords[:,0] + faceCoords[:,1]) / 2.
            >>> numerix.allclose(faceCenters, mesh.getFaceCenters(), atol = 1e-10, rtol = 1e-10)
            1

            >>> xc = dy  / Numeric.sqrt((dx * dx) + (dy * dy))
            >>> yc = dx  / Numeric.sqrt((dx * dx) + (dy * dy))
            
            >>> faceNormals = Numeric.array((( 0.0,-1.0), (0.0,-1.0), (0.0,-1.0),
            ...                              ( 0.0, 1.0), (0.0, 1.0), (0.0, 1.0),
            ...                              ( 0.0, 1.0), (0.0, 1.0), (0.0, 1.0),
            ...                              (-1.0, 0.0), (1.0, 0.0), (1.0, 0.0), (1.0, 0.0),
            ...                              (-1.0, 0.0), (1.0, 0.0), (1.0, 0.0), (1.0, 0.0),
            ...                              (  xc, -yc), ( xc, -yc), ( xc, -yc), ( xc, -yc), ( xc, -yc), ( xc, -yc),
            ...                              ( -xc, -yc), (-xc, -yc), (-xc, -yc), (-xc, -yc), (-xc, -yc), (-xc, -yc),
            ...                              ( -xc, -yc), (-xc, -yc), (-xc, -yc), (-xc, -yc), (-xc, -yc), (-xc, -yc),
            ...                              ( -xc,  yc), (-xc,  yc), (-xc,  yc), (-xc,  yc), (-xc,  yc), (-xc,  yc)))
            >>> numerix.allclose(faceNormals, mesh._getFaceNormals(), atol = 1e-10, rtol = 1e-10)
            1

            >>> cellToFaceOrientations = Numeric.array(((1, 1, 1), ( 1, 1, 1), ( 1, 1, 1), ( 1, 1, 1), ( 1, 1, 1), ( 1, 1, 1),
            ...                                         (1, 1,-1), ( 1, 1,-1), ( 1, 1,-1), ( 1, 1,-1), ( 1, 1,-1), ( 1, 1,-1),
            ...                                         (1, 1,-1), (-1, 1,-1), (-1, 1,-1), ( 1, 1,-1), (-1, 1,-1), (-1, 1,-1),
            ...                                         (1,-1,-1), ( 1,-1,-1), ( 1,-1,-1), (-1,-1,-1), (-1,-1,-1), (-1,-1,-1)))
            >>> numerix.allequal(cellToFaceOrientations, mesh._getCellFaceOrientations())
            1
                                             
            >>> cellVolumes = Numeric.array((0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
            ...                              0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
            ...                              0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25))
            >>> numerix.allclose(cellVolumes, mesh.getCellVolumes(), atol = 1e-10, rtol = 1e-10)
            1

            >>> sixth = 1.0 / 6.0
            >>> cellCenters = Numeric.array(((5*sixth, 0.5), (11*sixth, 0.5), (17*sixth, 0.5), (5*sixth, 1.5), (11*sixth, 1.5), (17*sixth, 1.5),
            ...                              (0.5, 5*sixth), (1.5, 5*sixth), (2.5, 5*sixth), (0.5, 11*sixth), (1.5, 11*sixth), (2.5, 11*sixth),
            ...                              (1*sixth, 0.5), (7*sixth, 0.5), (13*sixth, 0.5), (1*sixth, 1.5), (7*sixth, 1.5), (13*sixth, 1.5),
            ...                              (0.5, 1*sixth), (1.5, 1*sixth), (2.5, 1*sixth), (0.5, 7*sixth), (1.5, 7*sixth), (2.5, 7*sixth)))
            >>> cellCenters *= Numeric.array((dx, dy))
            >>> numerix.allclose(cellCenters, mesh.getCellCenters(), atol = 1e-10, rtol = 1e-10)
            1
                                              
            >>> yd = Numeric.sqrt(((dx/12.0)*(dx/12.0)) + ((dy/ 4.0)*(dy/ 4.0)))
            >>> xd = Numeric.sqrt(((dx/ 4.0)*(dx/ 4.0)) + ((dy/12.0)*(dy/12.0)))
            >>> faceToCellDistances = MA.masked_values(((dy/6.0, -1), (dy/6.0, -1), (dy/6.0, -1), (dy/6.0, dy/6.0), (dy/6.0, dy/6.0), (dy/6.0, dy/6.0), (dy/6.0, -1), (dy/6.0, -1), (dy/6.0, -1),
            ...                                         (dx/6.0, -1), (dx/6.0, dx/6.0), (dx/6.0, dx/6.0), (dx/6.0, -1), (dx/6.0, -1), (dx/6.0, dx/6.0), (dx/6.0, dx/6.0), (dx/6.0, -1),
            ...                                         (yd, xd), (yd, xd), (yd, xd), (yd, xd), (yd, xd), (yd, xd),
            ...                                         (yd, xd), (yd, xd), (yd, xd), (yd, xd), (yd, xd), (yd, xd),
            ...                                         (xd, yd), (xd, yd), (xd, yd), (xd, yd), (xd, yd), (xd, yd),
            ...                                         (yd, xd), (yd, xd), (yd, xd), (yd, xd), (yd, xd), (yd, xd)), -1)
            >>> numerix.allclose(faceToCellDistances, mesh._getFaceToCellDistances(), atol = 1e-10, rtol = 1e-10)
            1
                                              
            >>> dd = Numeric.sqrt((dx*dx)+(dy*dy)) / 3.0
            >>> cellDistances = Numeric.array((dy/6.0, dy/6.0, dy/6.0, dy/3.0, dy/3.0, dy/3.0, dy/6.0, dy/6.0, dy/6.0,
            ...                                dx/6.0, dx/3.0, dx/3.0, dx/6.0, dx/6.0, dx/3.0, dx/3.0, dx/6.0,
            ...                                dd, dd, dd, dd, dd, dd, dd, dd, dd, dd, dd, dd,
            ...                                dd, dd, dd, dd, dd, dd, dd, dd, dd, dd, dd, dd))
            >>> numerix.allclose(cellDistances, mesh._getCellDistances(), atol = 1e-10, rtol = 1e-10)
            1
            
            >>> faceToCellDistanceRatios = faceToCellDistances[...,0] / cellDistances
            >>> numerix.allclose(faceToCellDistanceRatios, mesh._getFaceToCellDistanceRatio(), atol = 1e-10, rtol = 1e-10)
            1

            >>> areaProjections = faceNormals * faceAreas[...,Numeric.NewAxis]
            >>> numerix.allclose(areaProjections, mesh._getAreaProjections(), atol = 1e-10, rtol = 1e-10)
            1

            >>> tangents1 = Numeric.array(((1.0, 0.0), (1.0, 0.0), (1.0, 0.0),
            ...                            (-1.0, 0.0), (-1.0, 0.0), (-1.0, 0.0),
            ...                            (-1.0, 0.0), (-1.0, 0.0), (-1.0, 0.0),
            ...                            (0.0, -1.0), (0.0, 1.0), (0.0, 1.0), (0.0, 1.0),
            ...                            (0.0, -1.0), (0.0, 1.0), (0.0, 1.0), (0.0, 1.0),
            ...                            (yc, xc),(yc, xc),(yc, xc),(yc, xc),(yc, xc),(yc, xc),
            ...                            (yc, -xc),(yc, -xc),(yc, -xc),(yc, -xc),(yc, -xc),(yc, -xc),
            ...                            (yc, -xc),(yc, -xc),(yc, -xc),(yc, -xc),(yc, -xc),(yc, -xc),
            ...                            (-yc, -xc),(-yc, -xc),(-yc, -xc),(-yc, -xc),(-yc, -xc),(-yc, -xc)))
            >>> numerix.allclose(tangents1, mesh._getFaceTangents1(), atol = 1e-10, rtol = 1e-10)
            1

            >>> tangents2 = Numeric.zeros((41, 2))
            >>> numerix.allclose(tangents2, mesh._getFaceTangents2(), atol = 1e-10, rtol = 1e-10)
            1

            >>> cellToCellIDs = MA.masked_values(((13, 18, 6), (14, 19, 7), (-1, 20, 8), (16, 21, 9), (17, 22, 10), (-1, 23, 11),
            ...                                   (21, 12, 0), (22, 13, 1), (23, 14, 2), (-1, 15, 3), (-1, 16, 4), (-1, 17, 5),
            ...                                   (-1, 18, 6), (0, 19, 7), (1, 20, 8), (-1, 21, 9), (3, 22, 10), (4, 23, 11),
            ...                                   (-1, 12, 0), (-1, 13, 1), (-1, 14, 2), (6, 15, 3), (7, 16, 4), (8, 17, 5)), -1)
            >>> numerix.allequal(cellToCellIDs, mesh._getCellToCellIDs())
            1

            >>> cellToCellDistances = Numeric.array([[dx/3.0, dd, dd],
            ...                                      [dx/3.0, dd, dd],
            ...                                      [dx/6.0, dd, dd],
            ...                                      [dx/3.0, dd, dd],
            ...                                      [dx/3.0, dd, dd],
            ...                                      [dx/6.0, dd, dd],
            ...                                      [dy/3.0, dd, dd],
            ...                                      [dy/3.0, dd, dd],
            ...                                      [dy/3.0, dd, dd],
            ...                                      [dy/6.0, dd, dd],
            ...                                      [dy/6.0, dd, dd],
            ...                                      [dy/6.0, dd, dd],
            ...                                      [dx/6.0, dd, dd],
            ...                                      [dx/3.0, dd, dd],
            ...                                      [dx/3.0, dd, dd],
            ...                                      [dx/6.0, dd, dd],
            ...                                      [dx/3.0, dd, dd],
            ...                                      [dx/3.0, dd, dd],
            ...                                      [dy/6.0, dd, dd],
            ...                                      [dy/6.0, dd, dd],
            ...                                      [dy/6.0, dd, dd],
            ...                                      [dy/3.0, dd, dd],
            ...                                      [dy/3.0, dd, dd],
            ...                                      [dy/3.0, dd, dd]])
            >>> numerix.allclose(cellToCellDistances, mesh._getCellToCellDistances(), atol = 1e-10, rtol = 1e-10)
            1

            >>> interiorCellIDs = Numeric.array((0, 1, 3, 4, 6, 7, 8, 13, 14, 16, 17, 21, 22, 23))
            >>> numerix.allequal(interiorCellIDs, mesh._getInteriorCellIDs())
            1

            >>> exteriorCellIDs = Numeric.array((2, 5, 9, 10, 11, 12, 15, 18, 19, 20))
            >>> numerix.allequal(exteriorCellIDs, mesh._getExteriorCellIDs())
            1

            >>> cellNormals = Numeric.array((((1, 0), (-xc, -yc), (-xc, yc)),
            ...                              ((1, 0), (-xc, -yc), (-xc, yc)),
            ...                              ((1, 0), (-xc, -yc), (-xc, yc)),
            ...                              ((1, 0), (-xc, -yc), (-xc, yc)),
            ...                              ((1, 0), (-xc, -yc), (-xc, yc)),
            ...                              ((1, 0), (-xc, -yc), (-xc, yc)),
            ...                              ((0, 1), (-xc, -yc), (xc, -yc)),
            ...                              ((0, 1), (-xc, -yc), (xc, -yc)),
            ...                              ((0, 1), (-xc, -yc), (xc, -yc)),
            ...                              ((0, 1), (-xc, -yc), (xc, -yc)),
            ...                              ((0, 1), (-xc, -yc), (xc, -yc)),
            ...                              ((0, 1), (-xc, -yc), (xc, -yc)),
            ...                              ((-1, 0), (xc, -yc), (xc, yc)),
            ...                              ((-1, 0), (xc, -yc), (xc, yc)),
            ...                              ((-1, 0), (xc, -yc), (xc, yc)),
            ...                              ((-1, 0), (xc, -yc), (xc, yc)),
            ...                              ((-1, 0), (xc, -yc), (xc, yc)),
            ...                              ((-1, 0), (xc, -yc), (xc, yc)),
            ...                              ((0, -1), (-xc, yc), (xc, yc)),
            ...                              ((0, -1), (-xc, yc), (xc, yc)),
            ...                              ((0, -1), (-xc, yc), (xc, yc)),
            ...                              ((0, -1), (-xc, yc), (xc, yc)),
            ...                              ((0, -1), (-xc, yc), (xc, yc)),
            ...                              ((0, -1), (-xc, yc), (xc, yc)) ))
            >>> numerix.allclose(cellNormals, mesh._getCellNormals(), atol = 1e-10, rtol = 1e-10)
            1

            >>> cellAreaProjections = Numeric.array(((( dy,  0), (-dy/2.,-dx/2.), (-dy/2., dx/2.)),
            ...                                      (( dy,  0), (-dy/2.,-dx/2.), (-dy/2., dx/2.)),
            ...                                      (( dy,  0), (-dy/2.,-dx/2.), (-dy/2., dx/2.)),
            ...                                      (( dy,  0), (-dy/2.,-dx/2.), (-dy/2., dx/2.)),
            ...                                      (( dy,  0), (-dy/2.,-dx/2.), (-dy/2., dx/2.)),
            ...                                      (( dy,  0), (-dy/2.,-dx/2.), (-dy/2., dx/2.)),
            ...                                      ((  0, dx), (-dy/2.,-dx/2.), ( dy/2.,-dx/2.)),
            ...                                      ((  0, dx), (-dy/2.,-dx/2.), ( dy/2.,-dx/2.)),
            ...                                      ((  0, dx), (-dy/2.,-dx/2.), ( dy/2.,-dx/2.)),
            ...                                      ((  0, dx), (-dy/2.,-dx/2.), ( dy/2.,-dx/2.)),
            ...                                      ((  0, dx), (-dy/2.,-dx/2.), ( dy/2.,-dx/2.)),
            ...                                      ((  0, dx), (-dy/2.,-dx/2.), ( dy/2.,-dx/2.)),
            ...                                      ((-dy,  0), ( dy/2.,-dx/2.), ( dy/2., dx/2.)),
            ...                                      ((-dy,  0), ( dy/2.,-dx/2.), ( dy/2., dx/2.)),
            ...                                      ((-dy,  0), ( dy/2.,-dx/2.), ( dy/2., dx/2.)),
            ...                                      ((-dy,  0), ( dy/2.,-dx/2.), ( dy/2., dx/2.)),
            ...                                      ((-dy,  0), ( dy/2.,-dx/2.), ( dy/2., dx/2.)),
            ...                                      ((-dy,  0), ( dy/2.,-dx/2.), ( dy/2., dx/2.)),
            ...                                      ((  0,-dx), (-dy/2., dx/2.), ( dy/2., dx/2.)),
            ...                                      ((  0,-dx), (-dy/2., dx/2.), ( dy/2., dx/2.)),
            ...                                      ((  0,-dx), (-dy/2., dx/2.), ( dy/2., dx/2.)),
            ...                                      ((  0,-dx), (-dy/2., dx/2.), ( dy/2., dx/2.)),
            ...                                      ((  0,-dx), (-dy/2., dx/2.), ( dy/2., dx/2.)),
            ...                                      ((  0,-dx), (-dy/2., dx/2.), ( dy/2., dx/2.)) ))
            >>> numerix.allclose(cellAreaProjections, mesh._getCellAreaProjections(), atol = 1e-10, rtol = 1e-10)
            1

            >>> tmp1 = Numeric.array((12, 5, 1))
            >>> tmp2 = Numeric.array((12, 5, 4))
            >>> tmp3 = Numeric.array((12, 4, 0))
            >>> tmp4 = Numeric.array((12, 1, 0))
            >>> tmp5 = Numeric.array((0, 1, 1))
            >>> cellVertexIDs = Numeric.array((tmp1, tmp1 + 1, tmp1 + 2, tmp1 + 3 + tmp5, tmp1 + 4 + tmp5, tmp1 + 5 + tmp5,
            ...                                tmp2, tmp2 + 1, tmp2 + 2, tmp2 + 3 + tmp5, tmp2 + 4 + tmp5, tmp2 + 5 + tmp5,
            ...                                tmp3, tmp3 + 1, tmp3 + 2, tmp3 + 3 + tmp5, tmp3 + 4 + tmp5, tmp3 + 5 + tmp5,
            ...                                tmp4, tmp4 + 1, tmp4 + 2, tmp4 + 3 + tmp5, tmp4 + 4 + tmp5, tmp4 + 5 + tmp5))

            >>> numerix.allclose(mesh._getCellVertexIDs(), cellVertexIDs)
            1
            

            >>> from fipy.tools import dump            
            >>> (f, filename) = dump.write(mesh, extension = '.gz')            
            >>> unpickledMesh = dump.read(filename, f)

            >>> numerix.allequal(mesh.getCellCenters(), unpickledMesh.getCellCenters())
            1
        """

## test test test
        
def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__": 
    _test() 







