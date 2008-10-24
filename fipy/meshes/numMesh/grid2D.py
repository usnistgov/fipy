#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "grid2D.py"
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
2D rectangular Mesh
"""
__docformat__ = 'restructuredtext'


from fipy.tools import numerix

from fipy.tools.inline import inline
from fipy.meshes.numMesh.mesh2D import Mesh2D
from fipy.tools import vector
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
        return "%s(dx=%s, dy=%s, nx=%d, ny=%d)" \
            % (self.__class__.__name__, `self.dx`, `self.dy`, self.nx, self.ny)
            
    def _createVertices(self):
        x = self._calcVertexCoordinates(self.dx, self.nx)
        x = numerix.resize(x, (self.numberOfVertices,))
            
        y = self._calcVertexCoordinates(self.dy, self.ny)
        y = numerix.repeat(y, self.nx + 1)
        
        return numerix.array((x, y))
    
    def _createFaces(self):
        """
        v1, v2 refer to the vertices.
        Horizontal faces are first
        """
        v1 = numerix.arange(self.numberOfVertices)
        v2 = v1 + 1

        horizontalFaces = vector.prune(numerix.array((v1, v2)), self.nx + 1, self.nx, axis=1)

        v1 = numerix.arange(self.numberOfVertices - (self.nx + 1))
        v2 = v1 + self.nx + 1
        verticalFaces =  numerix.array((v1, v2))

        ## The cell normals must point out of the cell.
        ## The left and bottom faces have only one neighboring cell,
        ## in the 2nd neighbor position (there is nothing in the 1st).
        ## 
        ## reverse some of the face orientations to obtain the correct normals

        tmp = horizontalFaces.copy()
        horizontalFaces[0,:self.nx] = tmp[1,:self.nx]
        horizontalFaces[1,:self.nx] = tmp[0,:self.nx]

        tmp = verticalFaces.copy()
        verticalFaces[0, :] = tmp[1, :]
        verticalFaces[1, :] = tmp[0, :]
        verticalFaces[0, ::(self.nx + 1)] = tmp[0, ::(self.nx + 1)]
        verticalFaces[1, ::(self.nx + 1)] = tmp[1,::(self.nx + 1)]

        return numerix.concatenate((horizontalFaces, verticalFaces), axis=1)

    def _createCells(self):
        """
        cells = (f1, f2, f3, f4) going anticlock wise.
        f1 etc. refer to the faces
        """
        self.numberOfHorizontalFaces = self.nx * (self.ny + 1)
        self.numberOfFaces = self.numberOfHorizontalFaces + self.ny * (self.nx + 1)
        return inline._optionalInline(self._createCellsIn, self._createCellsPy)

    def _createCellsPy(self):
        cellFaceIDs = numerix.zeros((4, self.nx * self.ny))
        faceIDs = numerix.arange(self.numberOfFaces)
        cellFaceIDs[0,:] = faceIDs[:self.numberOfHorizontalFaces - self.nx]
        cellFaceIDs[2,:] = cellFaceIDs[0,:] + self.nx
        cellFaceIDs[1,:] = vector.prune(faceIDs[self.numberOfHorizontalFaces:], self.nx + 1)
        cellFaceIDs[3,:] = cellFaceIDs[1,:] - 1
        return cellFaceIDs

    def _createCellsIn(self):
        cellFaceIDs = numerix.zeros((4, self.nx * self.ny))
        
        inline._runInline("""
            int ID = j * ni + i;
            int NCELLS = ni * nj;
            cellFaceIDs[ID + 0 * NCELLS] = ID;
            cellFaceIDs[ID + 2 * NCELLS] = cellFaceIDs[ID + 0 * NCELLS] + ni;
            cellFaceIDs[ID + 3 * NCELLS] = horizontalFaces + ID + j;
            cellFaceIDs[ID + 1 * NCELLS] = cellFaceIDs[ID + 3 * NCELLS] + 1;
        """,
        horizontalFaces=self.numberOfHorizontalFaces,
        cellFaceIDs=cellFaceIDs,
        ni=self.nx,
        nj=self.ny)

        return cellFaceIDs
    
    def getScale(self):
        return self.scale['length']
        
    def getPhysicalShape(self):
        """Return physical dimensions of Grid2D.
        """
        return PhysicalField(value = (self.nx * self.dx * self.getScale(), self.ny * self.dy * self.getScale()))

    def _getMeshSpacing(self):
        return numerix.array((self.dx,self.dy))[...,numerix.newaxis]
    
    def getShape(self):
        return (self.nx, self.ny)

    def _isOrthogonal(self):
        return True
    
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
            
            >>> vertices = numerix.array(((0., 1., 2., 3., 0., 1., 2., 3., 0., 1., 2., 3.),
            ...                           (0., 0., 0., 0., 1., 1., 1., 1., 2., 2., 2., 2.)))
            >>> vertices *= numerix.array(((dx,), (dy,)))
            >>> numerix.allequal(vertices, mesh._createVertices())
            1
        
            >>> faces = numerix.array(((1, 2, 3, 4, 5, 6, 8, 9, 10, 0, 5, 6, 7, 4, 9, 10, 11),
            ...                        (0, 1, 2, 5, 6, 7, 9, 10, 11, 4, 1, 2, 3, 8, 5, 6, 7)))
            >>> numerix.allequal(faces, mesh._createFaces())
            1

            >>> cells = numerix.array(((0, 1, 2, 3, 4, 5),
            ...                        (10, 11, 12, 14, 15, 16),
            ...                        (3, 4, 5, 6, 7, 8),
            ...                        (9, 10, 11, 13, 14, 15)))
            >>> numerix.allequal(cells, mesh._createCells())
            1

            >>> externalFaces = numerix.array((0, 1, 2, 6, 7, 8, 9 , 12, 13, 16))
            >>> numerix.allequal(externalFaces, 
            ...                  numerix.nonzero(mesh.getExteriorFaces()))
            1

            >>> internalFaces = numerix.array((3, 4, 5, 10, 11, 14, 15))
            >>> print numerix.allequal(internalFaces, 
            ...                        numerix.nonzero(mesh.getInteriorFaces()))
            1

            >>> from fipy.tools.numerix import MA
            >>> faceCellIds = MA.masked_values(((0, 1, 2, 0, 1, 2, 3, 4, 5, 0, 0, 1, 2, 3, 3, 4, 5),
            ...                                 (-1, -1, -1, 3, 4, 5, -1, -1, -1, -1, 1, 2, -1, -1, 4, 5, -1)), -1)
            >>> numerix.allequal(faceCellIds, mesh.getFaceCellIDs())
            1
            
            >>> faceAreas = numerix.array((dx, dx, dx, dx, dx, dx, dx, dx, dx,
            ...                            dy, dy, dy, dy, dy, dy, dy, dy))
            >>> numerix.allclose(faceAreas, mesh._getFaceAreas(), atol = 1e-10, rtol = 1e-10)
            1
            
            >>> faceCoords = numerix.take(vertices, faces, axis=1)
            >>> faceCenters = (faceCoords[...,0,:] + faceCoords[...,1,:]) / 2.
            >>> numerix.allclose(faceCenters, mesh.getFaceCenters(), atol = 1e-10, rtol = 1e-10)
            1

            >>> faceNormals = numerix.array(((0., 0., 0., 0., 0., 0., 0., 0., 0., -1., 1., 1., 1., -1., 1., 1., 1.),
            ...                              (-1., -1., -1., 1., 1., 1., 1., 1., 1., 0., 0., 0., 0., 0., 0., 0., 0.)))
            >>> numerix.allclose(faceNormals, mesh._getFaceNormals(), atol = 1e-10, rtol = 1e-10)
            1

            >>> cellToFaceOrientations = numerix.array(((1,  1,  1, -1, -1, -1),
            ...                                         (1,  1,  1,  1,  1,  1),
            ...                                         (1,  1,  1,  1,  1,  1),
            ...                                         (1, -1, -1,  1, -1, -1)))
            >>> numerix.allequal(cellToFaceOrientations, mesh._getCellFaceOrientations())
            1
                                             
            >>> cellVolumes = numerix.array((dx*dy, dx*dy, dx*dy, dx*dy, dx*dy, dx*dy))
            >>> numerix.allclose(cellVolumes, mesh.getCellVolumes(), atol = 1e-10, rtol = 1e-10)
            1

            >>> cellCenters = numerix.array(((dx/2., 3.*dx/2., 5.*dx/2., dx/2., 3.*dx/2., 5.*dx/2.),
            ...                              (dy/2., dy/2., dy/2., 3.*dy/2., 3.*dy/2., 3.*dy/2.)))
            >>> numerix.allclose(cellCenters, mesh.getCellCenters(), atol = 1e-10, rtol = 1e-10)
            1
                                              
            >>> faceToCellDistances = MA.masked_values(((dy / 2., dy / 2., dy / 2., dy / 2., dy / 2., dy / 2., dy / 2., dy / 2., dy / 2., dx / 2., dx / 2., dx / 2., dx / 2., dx / 2., dx / 2., dx / 2., dx / 2.),
            ...                                         (-1, -1, -1, dy / 2., dy / 2., dy / 2., -1, -1, -1, -1, dx / 2., dx / 2., -1, -1, dx / 2., dx / 2., -1)), -1)
            >>> numerix.allclose(faceToCellDistances, mesh._getFaceToCellDistances(), atol = 1e-10, rtol = 1e-10)
            1
                                              
            >>> cellDistances = numerix.array((dy / 2., dy / 2., dy / 2.,
            ...                                dy, dy, dy,
            ...                                dy / 2., dy / 2., dy / 2.,
            ...                                dx / 2., dx, dx,
            ...                                dx / 2.,
            ...                                dx / 2., dx, dx,
            ...                                dx / 2.))
            >>> numerix.allclose(cellDistances, mesh._getCellDistances(), atol = 1e-10, rtol = 1e-10)
            1
            
            >>> faceToCellDistanceRatios = faceToCellDistances[0] / cellDistances
            >>> numerix.allclose(faceToCellDistanceRatios, mesh._getFaceToCellDistanceRatio(), atol = 1e-10, rtol = 1e-10)
            1

            >>> areaProjections = faceNormals * faceAreas
            >>> numerix.allclose(areaProjections, mesh._getAreaProjections(), atol = 1e-10, rtol = 1e-10)
            1

            >>> tangents1 = numerix.array(((1., 1., 1., -1., -1., -1., -1., -1., -1., 0., 0., 0., 0., 0., 0., 0., 0.),
            ...                            (0., 0., 0., 0., 0., 0., 0., 0., 0., -1., 1., 1., 1., -1., 1., 1., 1.)))
            >>> numerix.allclose(tangents1, mesh._getFaceTangents1(), atol = 1e-10, rtol = 1e-10)
            1

            >>> tangents2 = numerix.zeros((2, 17), 'd')
            >>> numerix.allclose(tangents2, mesh._getFaceTangents2(), atol = 1e-10, rtol = 1e-10)
            1

            >>> cellToCellIDs = MA.masked_values(((-1, -1, -1, 0, 1, 2),
            ...                                   (1, 2, -1, 4, 5, -1),
            ...                                   (3, 4, 5, -1, -1, -1),
            ...                                   (-1, 0, 1, -1, 3, 4)), -1)
            >>> numerix.allequal(cellToCellIDs, mesh._getCellToCellIDs())
            1

            >>> cellToCellDistances = MA.masked_values(((dy / 2., dy / 2., dy / 2.,      dy,      dy,      dy),
            ...                                         (     dx,      dx, dx / 2.,      dx,      dx, dx / 2.),
            ...                                         (     dy,      dy,      dy, dy / 2., dy / 2., dy / 2.),
            ...                                         (dx / 2.,      dx,      dx, dx / 2.,      dx,      dx)), -1)
            >>> numerix.allclose(cellToCellDistances, mesh._getCellToCellDistances(), atol = 1e-10, rtol = 1e-10)
            1

            >>> interiorCellIDs = numerix.array(())
            >>> numerix.allequal(interiorCellIDs, mesh._getInteriorCellIDs())
            1

            >>> exteriorCellIDs = numerix.array((0, 1, 2, 3, 4, 5))
            >>> numerix.allequal(exteriorCellIDs, mesh._getExteriorCellIDs())
            1

            >>> cellNormals = numerix.array(((( 0,  0,  0,  0,  0,  0),
            ...                               ( 1,  1,  1,  1,  1,  1),
            ...                               ( 0,  0,  0,  0,  0,  0),
            ...                               (-1, -1, -1, -1, -1, -1)),
            ...                              ((-1, -1, -1, -1, -1, -1),
            ...                               ( 0,  0,  0,  0,  0,  0),
            ...                               ( 1,  1,  1,  1,  1,  1),
            ...                               ( 0,  0,  0,  0,  0,  0))))
            >>> numerix.allclose(cellNormals, mesh._getCellNormals(), atol = 1e-10, rtol = 1e-10)
            1

            >>> cellAreaProjections = numerix.array((((  0,  0,  0,  0,  0,  0),
            ...                                       ( dy, dy, dy, dy, dy, dy),
            ...                                       (  0,  0,  0,  0,  0,  0),
            ...                                       (-dy,-dy,-dy,-dy,-dy,-dy)),
            ...                                      ((-dx,-dx,-dx,-dx,-dx,-dx),
            ...                                       (  0,  0,  0,  0,  0,  0),
            ...                                       ( dx, dx, dx, dx, dx, dx),
            ...                                       (  0,  0,  0,  0,  0,  0))))
            >>> numerix.allclose(cellAreaProjections, mesh._getCellAreaProjections(), atol = 1e-10, rtol = 1e-10)
            1

            >>> cellVertexIDs = MA.masked_values(((5, 6, 7, 9, 10, 11),
            ...                                   (4, 5, 6, 8,  9, 10),
            ...                                   (1, 2, 3, 5,  6,  7),
            ...                                   (0, 1, 2, 4,  5,  6)), -1000)

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
