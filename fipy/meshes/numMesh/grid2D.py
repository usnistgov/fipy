#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "grid2D.py"
 #                                    created: 11/10/03 {3:30:42 PM} 
 #                                last update: 5/15/06 {4:17:38 PM} 
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
2D rectangular Mesh
"""
__docformat__ = 'restructuredtext'


import Numeric

from fipy.tools.inline import inline
from fipy.meshes.numMesh.mesh2D import Mesh2D
from fipy.meshes.meshIterator import FaceIterator
from fipy.tools import vector
from fipy.tools import numerix
from fipy.tools.dimensions.physicalField import PhysicalField

class Grid2D(Mesh2D):
    """
    Creates a 2D grid mesh with horizontal faces numbered
    first and then vertical faces.
    """
    def __init__(self, dx = 1., dy = 1., nx = None, ny = None):
	self.dx = PhysicalField(value = dx)
	scale = PhysicalField(value = 1, unit = self.dx.getUnit())
	self.dx /= scale
        
        self.nx = self._calcNumPts(d = self.dx, n = nx, axis = "x")
        
	self.dy = PhysicalField(value = dy)
	if self.dy.getUnit().isDimensionless():
	    self.dy = dy
	else:
	    self.dy /= scale
            
        self.ny = self._calcNumPts(d = self.dy, n = ny, axis = "y")
	
        self.numberOfVertices = (self.nx + 1) * (self. ny + 1)
	
	vertices = self._createVertices()
        faces = self._createFaces()
        cells = self._createCells()
        Mesh2D.__init__(self, vertices, faces, cells)
	
	self.setScale(value = scale)
        
    def __repr__(self):
        return "%s(dx = %s, dy = %s, nx = %d, ny = %d)" \
            % (self.__class__.__name__, `self.dx`, `self.dy`, self.nx, self.ny)
            
    def _createVertices(self):
        x = self._calcVertexCoordinates(self.dx, self.nx)
        x = Numeric.resize(x, (self.numberOfVertices,))
            
        y = self._calcVertexCoordinates(self.dy, self.ny)
        y = Numeric.repeat(y, self.nx + 1)
        
        return Numeric.transpose(Numeric.array((x, y)))
    
    def _createFaces(self):
        """
        v1, v2 refer to the vertices.
        Horizontal faces are first
        """
        v1 = Numeric.arange(self.numberOfVertices)
        v2 = v1 + 1
        horizontalFaces = vector.prune(Numeric.transpose(Numeric.array((v1, v2))), self.nx + 1, self.nx)
        v1 = Numeric.arange(self.numberOfVertices - (self.nx + 1))
        v2 = v1 + self.nx + 1
        verticalFaces =  Numeric.transpose(Numeric.array((v1, v2)))

        ## The cell normals must point out of the cell.
        ## The left and bottom faces have only one neighboring cell,
        ## in the 2nd neighbor position (there is nothing in the 1st).
        ## 
        ## reverse some of the face orientations to obtain the correct normals

        tmp = horizontalFaces.copy()
        horizontalFaces[:self.nx, 0] = tmp[:self.nx, 1]
        horizontalFaces[:self.nx, 1] = tmp[:self.nx, 0]

        tmp = verticalFaces.copy()
        verticalFaces[:, 0] = tmp[:, 1]
        verticalFaces[:, 1] = tmp[:, 0]
        verticalFaces[::(self.nx + 1), 0] = tmp[::(self.nx + 1), 0]
        verticalFaces[::(self.nx + 1), 1] = tmp[::(self.nx + 1), 1]

        return Numeric.concatenate((horizontalFaces, verticalFaces))

    def _createCells(self):
        """
        cells = (f1, f2, f3, f4) going anticlock wise.
        f1 etc. refer to the faces
        """
        self.numberOfHorizontalFaces = self.nx * (self.ny + 1)
        self.numberOfFaces = self.numberOfHorizontalFaces + self.ny * (self.nx + 1)
        return inline._optionalInline(self._createCellsIn, self._createCellsPy)

    def _createCellsPy(self):
        cellFaceIDs = Numeric.zeros((self.nx * self.ny, 4))
        faceIDs = Numeric.arange(self.numberOfFaces)
        cellFaceIDs[:,0] = faceIDs[:self.numberOfHorizontalFaces - self.nx]
        cellFaceIDs[:,2] = cellFaceIDs[:,0] + self.nx
        cellFaceIDs[:,1] = vector.prune(faceIDs[self.numberOfHorizontalFaces:], self.nx + 1)
        cellFaceIDs[:,3] = cellFaceIDs[:,1] - 1
        return cellFaceIDs

    def _createCellsIn(self):
        cellFaceIDs = Numeric.zeros((self.nx * self.ny, 4))
        
        inline._runInlineLoop2("""
            int ID = j * ni + i;
            cellFaceIDs(ID, 0) = ID;
            cellFaceIDs(ID, 2) = cellFaceIDs(ID, 0) + ni;
            cellFaceIDs(ID, 3) = horizontalFaces + ID + j;
            cellFaceIDs(ID, 1) = cellFaceIDs(ID, 3) + 1;
	""",
        horizontalFaces=self.numberOfHorizontalFaces,
        cellFaceIDs=cellFaceIDs,
        ni=self.nx,
        nj=self.ny)

        return cellFaceIDs

    def getFacesLeft(self):
	"""
        Return list of faces on left boundary of Grid2D.
        
            >>> mesh = Grid2D(nx = 3, ny = 2, dx = 0.5, dy = 2.)        
            >>> numerix.allequal((9, 13), mesh.getFacesLeft())
            1
	"""
        return FaceIterator(mesh = self, ids = Numeric.arange(self.numberOfHorizontalFaces, self.numberOfFaces, self.nx + 1))
	
    def getFacesRight(self):
	"""
        Return list of faces on right boundary of Grid2D.
        
            >>> mesh = Grid2D(nx = 3, ny = 2, dx = 0.5, dy = 2.)        
            >>> numerix.allequal((12, 16), mesh.getFacesRight())
            1
	"""
        return FaceIterator(mesh = self, ids = Numeric.arange(self.numberOfHorizontalFaces + self.nx, self.numberOfFaces, self.nx + 1))
	
    def getFacesTop(self):
	"""
        Return list of faces on top boundary of Grid2D.
        
            >>> mesh = Grid2D(nx = 3, ny = 2, dx = 0.5, dy = 2.)        
            >>> numerix.allequal((6, 7, 8), mesh.getFacesTop())
            1
	"""
        return FaceIterator(mesh = self, ids = Numeric.arange(self.numberOfHorizontalFaces - self.nx, self.numberOfHorizontalFaces))
	
    def getFacesBottom(self):
	"""
        Return list of faces on bottom boundary of Grid2D.
        
            >>> mesh = Grid2D(nx = 3, ny = 2, dx = 0.5, dy = 2.)        
            >>> numerix.allequal((0, 1, 2), mesh.getFacesBottom())
            1
	"""
        return FaceIterator(mesh = self, ids = Numeric.arange(self.nx))
        
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
        """
        Used internally to collect the necessary information to ``pickle`` the 
        `Grid2D` to persistent storage.
        """
        return {
            'dx' : self.dx,            
            'dy' : self.dy,
            'nx' : self.nx,
            'ny' : self.ny
        }

    def __setstate__(self, dict):
        """
        Used internally to create a new `Grid2D` from ``pickled`` 
        persistent storage.
        """
        self.__init__(dx = dict['dx'], dy = dict['dy'], nx = dict['nx'], ny = dict['ny'])

    def _test(self):
        """
        These tests are not useful as documentation, but are here to ensure
        everything works as expected.

            >>> dx = 0.5
            >>> dy = 2.
            >>> nx = 3
            >>> ny = 2
            
            >>> mesh = Grid2D(nx = nx, ny = ny, dx = dx, dy = dy)     
            
            >>> vertices = Numeric.array(((0., 0.), (1., 0.), (2., 0.), (3., 0.),
            ...                           (0., 1.), (1., 1.), (2., 1.), (3., 1.),
            ...                           (0., 2.), (1., 2.), (2., 2.), (3., 2.)))
            >>> vertices *= Numeric.array((dx, dy))
            >>> numerix.allequal(vertices, mesh._createVertices())
            1
        
            >>> faces = Numeric.array(((1, 0), (2, 1), (3, 2),
            ...                        (4, 5), (5, 6), (6, 7),
            ...                        (8, 9), (9, 10), (10, 11),
            ...                        (0, 4), (5, 1), (6, 2), (7, 3),
            ...                        (4, 8), (9, 5), (10, 6), (11, 7)))
            >>> numerix.allequal(faces, mesh._createFaces())
            1

            >>> cells = Numeric.array(((0, 10, 3, 9),
            ...                       (1 , 11, 4, 10),
            ...                       (2, 12, 5, 11),
            ...                       (3, 14, 6, 13),
            ...                       (4, 15, 7, 14),
            ...                       (5, 16, 8, 15)))
            >>> numerix.allequal(cells, mesh._createCells())
            1

            >>> externalFaces = Numeric.array((0, 1, 2, 6, 7, 8, 9 , 12, 13, 16))
            >>> tmp = list(mesh.getExteriorFaces())
            >>> tmp.sort()
            >>> numerix.allequal(externalFaces, tmp)
            1

            >>> internalFaces = Numeric.array((3, 4, 5, 10, 11, 14, 15))
            >>> numerix.allequal(internalFaces, mesh.getInteriorFaces())
            1

            >>> import MA
            >>> faceCellIds = MA.masked_values(((0, -1), (1, -1), (2, -1),
            ...                                 (0, 3), (1, 4), (2, 5),
            ...                                 (3, -1), (4, -1), (5, -1),
            ...                                 (0, -1), (0, 1), (1, 2), (2, -1),
            ...                                 (3, -1), (3, 4), (4, 5), (5, -1)), -1)
            >>> numerix.allequal(faceCellIds, mesh.getFaceCellIDs())
            1
            
            >>> faceAreas = Numeric.array((dx, dx, dx, dx, dx, dx, dx, dx, dx,
            ...                            dy, dy, dy, dy, dy, dy, dy, dy))
            >>> numerix.allclose(faceAreas, mesh._getFaceAreas(), atol = 1e-10, rtol = 1e-10)
            1
            
            >>> faceCoords = Numeric.take(vertices, faces)
            >>> faceCenters = (faceCoords[:,0] + faceCoords[:,1]) / 2.
            >>> numerix.allclose(faceCenters, mesh.getFaceCenters(), atol = 1e-10, rtol = 1e-10)
            1

            >>> faceNormals = Numeric.array(((0., -1.), (0., -1.), (0., -1.),
            ...                              (0., 1.), (0., 1.), (0., 1.),
            ...                              (0., 1.), (0., 1.), (0., 1.),
            ...                              (-1., 0), (1., 0), (1., 0), (1., 0),
            ...                              (-1., 0), (1., 0), (1., 0), (1., 0)))
            >>> numerix.allclose(faceNormals, mesh._getFaceNormals(), atol = 1e-10, rtol = 1e-10)
            1

            >>> cellToFaceOrientations = Numeric.array(((1, 1, 1, 1), (1, 1, 1, -1), (1, 1, 1, -1),
            ...                                         (-1, 1, 1, 1), (-1, 1, 1, -1), (-1, 1, 1, -1)))
            >>> numerix.allequal(cellToFaceOrientations, mesh._getCellFaceOrientations())
            1
                                             
            >>> cellVolumes = Numeric.array((dx*dy, dx*dy, dx*dy, dx*dy, dx*dy, dx*dy))
            >>> numerix.allclose(cellVolumes, mesh.getCellVolumes(), atol = 1e-10, rtol = 1e-10)
            1

            >>> cellCenters = Numeric.array(((dx/2.,dy/2.), (3.*dx/2.,dy/2.), (5.*dx/2.,dy/2.),
            ...                              (dx/2.,3.*dy/2.), (3.*dx/2.,3.*dy/2.), (5.*dx/2.,3.*dy/2.)))
            >>> numerix.allclose(cellCenters, mesh.getCellCenters(), atol = 1e-10, rtol = 1e-10)
            1
                                              
            >>> faceToCellDistances = MA.masked_values(((dy / 2., -1), (dy / 2., -1), (dy / 2., -1),
            ...                                         (dy / 2., dy / 2.), (dy / 2., dy / 2.), (dy / 2., dy / 2.),
            ...                                         (dy / 2., -1), (dy / 2., -1), (dy / 2., -1),
            ...                                         (dx / 2., -1), (dx / 2., dx / 2.), (dx / 2., dx / 2.),
            ...                                         (dx / 2., -1),
            ...                                         (dx / 2., -1), (dx / 2., dx / 2.), (dx / 2., dx / 2.),
            ...                                         (dx / 2., -1)), -1)
            >>> numerix.allclose(faceToCellDistances, mesh._getFaceToCellDistances(), atol = 1e-10, rtol = 1e-10)
            1
                                              
            >>> cellDistances = Numeric.array((dy / 2., dy / 2., dy / 2.,
            ...                                dy, dy, dy,
            ...                                dy / 2., dy / 2., dy / 2.,
            ...                                dx / 2., dx, dx,
            ...                                dx / 2.,
            ...                                dx / 2., dx, dx,
            ...                                dx / 2.))
            >>> numerix.allclose(cellDistances, mesh._getCellDistances(), atol = 1e-10, rtol = 1e-10)
            1
            
            >>> faceToCellDistanceRatios = faceToCellDistances[...,0] / cellDistances
            >>> numerix.allclose(faceToCellDistanceRatios, mesh._getFaceToCellDistanceRatio(), atol = 1e-10, rtol = 1e-10)
            1

            >>> areaProjections = faceNormals * faceAreas[...,Numeric.NewAxis]
            >>> numerix.allclose(areaProjections, mesh._getAreaProjections(), atol = 1e-10, rtol = 1e-10)
            1

            >>> tangents1 = Numeric.array(((1., 0), (1., 0),(1., 0),
            ...                            (-1., 0), (-1., 0),(-1., 0),
            ...                            (-1., 0), (-1., 0),(-1., 0),
            ...                            (0., -1.), (0., 1.), (0., 1.), (0., 1.),
            ...                            (0., -1.), (0., 1.), (0., 1.), (0., 1.)))
            >>> numerix.allclose(tangents1, mesh._getFaceTangents1(), atol = 1e-10, rtol = 1e-10)
            1

            >>> tangents2 = Numeric.array(((0., 0), (0., 0),(0., 0),
            ...                            (-0., 0), (-0., 0),(-0., 0),
            ...                            (-0., 0), (-0., 0),(-0., 0),
            ...                            (0., -0.), (0., 0.), (0., 0.), (0., 0.),
            ...                            (0., -0.), (0., 0.), (0., 0.), (0., 0.)))
            >>> numerix.allclose(tangents2, mesh._getFaceTangents2(), atol = 1e-10, rtol = 1e-10)
            1

            >>> cellToCellIDs = MA.masked_values(((-1, 1, 3, -1),
            ...                                   (-1, 2, 4, 0),
            ...                                   (-1, -1, 5, 1),
            ...                                   (0, 4, -1, -1),
            ...                                   (1, 5, -1, 3),
            ...                                   (2, -1, -1, 4)), -1)
            >>> numerix.allequal(cellToCellIDs, mesh._getCellToCellIDs())
            1

            >>> cellToCellDistances = MA.masked_values(((dy / 2., dx, dy, dx / 2.),
            ...                                         (dy/ 2., dx, dy, dx),
            ...                                         (dy / 2., dx / 2., dy, dx),
            ...                                         (dy, dx, dy / 2., dx / 2.),
            ...                                         (dy, dx, dy / 2., dx),
            ...                                         (dy, dx / 2., dy / 2., dx)), -1)
            >>> numerix.allclose(cellToCellDistances, mesh._getCellToCellDistances(), atol = 1e-10, rtol = 1e-10)
            1

            >>> interiorCellIDs = Numeric.array(())
            >>> numerix.allequal(interiorCellIDs, mesh._getInteriorCellIDs())
            1

            >>> exteriorCellIDs = Numeric.array((0, 1, 2, 3, 4, 5))
            >>> numerix.allequal(exteriorCellIDs, mesh._getExteriorCellIDs())
            1

            >>> cellNormals = Numeric.array((((0, -1), (1, 0), (0, 1), (-1, 0)),
            ...                              ((0, -1), (1, 0), (0, 1), (-1, 0)),
            ...                              ((0, -1), (1, 0), (0, 1), (-1, 0)),
            ...                              ((0, -1), (1, 0), (0, 1), (-1, 0)),
            ...                              ((0, -1), (1, 0), (0, 1), (-1, 0)),
            ...                              ((0, -1), (1, 0), (0, 1), (-1, 0)) ))
            >>> numerix.allclose(cellNormals, mesh._getCellNormals(), atol = 1e-10, rtol = 1e-10)
            1

            >>> vv = Numeric.array(((0, -dx), (dy, 0), (0, dx), (-dy, 0)))
            >>> cellAreaProjections = Numeric.array(((vv,vv,vv,vv,vv,vv)))
            >>> numerix.allclose(cellAreaProjections, mesh._getCellAreaProjections(), atol = 1e-10, rtol = 1e-10)
            1

            >>> cellVertexIDs = MA.masked_array(((5, 4, 1, 0),
            ...                                  (6, 5, 2, 1),
            ...                                  (7, 6, 3, 2),
            ...                                  (9, 8, 5, 4),
            ...                                  (10, 9, 6, 5),
            ...                                  (11, 10, 7, 6)), -1000) 

            >>> numerix.allclose(mesh._getCellVertexIDs(), cellVertexIDs)
            1

            >>> from fipy.tools import dump            
            >>> (f, filename) = dump.write(mesh, extension = '.gz')            
            >>> unpickledMesh = dump.read(filename, f)

            >>> numerix.allequal(mesh.getCellCenters(), unpickledMesh.getCellCenters())
            1
        """

def _test():
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
