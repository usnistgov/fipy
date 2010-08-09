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

from fipy.meshes.numMesh.mesh2D import Mesh2D
from fipy.tools import inline
from fipy.tools import numerix
from fipy.tools import vector
from fipy.tools.dimensions.physicalField import PhysicalField
from fipy.tools import parallel

class Grid2D(Mesh2D):
    """
    Creates a 2D grid mesh with horizontal faces numbered
    first and then vertical faces.
    """
    def __init__(self, dx=1., dy=1., nx=None, ny=None, overlap=2, communicator=parallel):
        
        self.args = {
            'dx': dx, 
            'dy': dy, 
            'nx': nx, 
            'ny': ny, 
            'overlap': overlap,
            'communicator': communicator
        }

        self.dx = PhysicalField(value = dx)
        scale = PhysicalField(value=1, unit=self.dx.getUnit())
        self.dx /= scale
        
        nx = self._calcNumPts(d=self.dx, n=nx, axis="x")
        
        self.dy = PhysicalField(value = dy)
        if self.dy.getUnit().isDimensionless():
            self.dy = dy
        else:
            self.dy /= scale
            
        ny = self._calcNumPts(d=self.dy, n=ny, axis="y")
        
        (self.nx,
         self.ny,
         self.overlap,
         self.offset) = self._calcParallelGridInfo(nx, ny, overlap, communicator)

        if numerix.getShape(self.dx) is not ():
            Xoffset = numerix.sum(self.dx[0:self.offset[0]])
            self.dx = self.dx[self.offset[0]:self.offset[0] + self.nx]
        else:
            Xoffset = 0

        if numerix.getShape(self.dy) is not ():
            Yoffset =  numerix.sum(self.dy[0:self.offset[1]])
            self.dy = self.dy[self.offset[1]:self.offset[1] + self.ny]
        else:
            Yoffset = 0
            
        if self.nx == 0:
            self.ny = 0
        if self.ny == 0:
            self.nx = 0
        if self.nx == 0 or self.ny == 0:
            self.numberOfHorizontalRows = 0
            self.numberOfVerticalColumns = 0
        else:
            self.numberOfHorizontalRows = (self.ny + 1)
            self.numberOfVerticalColumns = (self.nx + 1)

        self.numberOfVertices = self.numberOfHorizontalRows * self.numberOfVerticalColumns

        vertices = self._createVertices() + ((Xoffset,), (Yoffset,))
        faces = self._createFaces()
        self.numberOfFaces = len(faces[0])
        cells = self._createCells()
        Mesh2D.__init__(self, vertices, faces, cells, communicator=communicator)
        
        self.setScale(value = scale)

    def _calcParallelGridInfo(self, nx, ny, overlap, communicator):
        
        procID = communicator.procID
        Nproc = communicator.Nproc

        overlap = min(overlap, ny)
        cellsPerNode = max(int(ny / Nproc), overlap)
        occupiedNodes = min(int(ny / (cellsPerNode or 1)), Nproc)
            
        overlap = {
            'left': 0,
            'right': 0,
            'bottom': overlap * (procID > 0) * (procID < occupiedNodes),
            'top': overlap * (procID < occupiedNodes - 1)
        }
        
        offset = (0,
                  min(procID, occupiedNodes-1) * cellsPerNode - overlap['bottom'])
                
        local_nx = nx
        local_ny = cellsPerNode * (procID < occupiedNodes)
        if procID == occupiedNodes - 1:
            local_ny += (ny - cellsPerNode * occupiedNodes)
        local_ny = local_ny + overlap['bottom'] + overlap['top']
        
        self.globalNumberOfCells = nx * ny
        self.globalNumberOfFaces = nx * (ny + 1) + ny * (nx + 1)
        
        return local_nx, local_ny, overlap, offset

    def __repr__(self):
        return "%s(dx=%s, dy=%s, nx=%d, ny=%d)" \
            % (self.__class__.__name__, str(self.args["dx"]), str(self.args["dy"]), 
               self.args["nx"], self.args["ny"])
            
    def _createVertices(self):
        x = self._calcVertexCoordinates(self.dx, self.nx)
        x = numerix.resize(x, (self.numberOfVertices,))
            
        y = self._calcVertexCoordinates(self.dy, self.ny)
        y = numerix.repeat(y, self.numberOfVerticalColumns)
        
        return numerix.array((x, y))
    
    def _createFaces(self):
        """
        v1, v2 refer to the vertices.
        Horizontal faces are first
        """
        v1 = numerix.arange(self.numberOfVertices)
        v2 = v1 + 1

        horizontalFaces = vector.prune(numerix.array((v1, v2)), self.numberOfVerticalColumns, self.nx, axis=1)

        v1 = numerix.arange(self.numberOfVertices - self.numberOfVerticalColumns)
        v2 = v1 + self.numberOfVerticalColumns
        verticalFaces =  numerix.array((v1, v2))

        ## The cell normals must point out of the cell.
        ## The left and bottom faces have only one neighboring cell,
        ## in the 2nd neighbor position (there is nothing in the 1st).
        ## 
        ## reverse some of the face orientations to obtain the correct normals

        tmp = horizontalFaces.copy()
        horizontalFaces[0,:self.nx] = tmp[1,:self.nx]
        horizontalFaces[1,:self.nx] = tmp[0,:self.nx]

        self.numberOfHorizontalFaces = horizontalFaces.shape[-1]

        tmp = verticalFaces.copy()
        verticalFaces[0, :] = tmp[1, :]
        verticalFaces[1, :] = tmp[0, :]
        if self.numberOfVerticalColumns > 0:
            verticalFaces[0, ::self.numberOfVerticalColumns] = tmp[0, ::self.numberOfVerticalColumns]
            verticalFaces[1, ::self.numberOfVerticalColumns] = tmp[1,::self.numberOfVerticalColumns]

        return numerix.concatenate((horizontalFaces, verticalFaces), axis=1)

    def _createCells(self):
        """
        cells = (f1, f2, f3, f4) going anticlock wise.
        f1 etc. refer to the faces
        """
        return inline._optionalInline(self._createCellsIn, self._createCellsPy)

    def _createCellsPy(self):
        cellFaceIDs = numerix.zeros((4, self.nx * self.ny))
        faceIDs = numerix.arange(self.numberOfFaces)
        if self.numberOfFaces > 0:
            cellFaceIDs[0,:] = faceIDs[:self.numberOfHorizontalFaces - self.nx]
            cellFaceIDs[2,:] = cellFaceIDs[0,:] + self.nx
            cellFaceIDs[1,:] = vector.prune(faceIDs[self.numberOfHorizontalFaces:], self.numberOfVerticalColumns)
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
        
    def _getGlobalNonOverlappingCellIDs(self):
        """
        Return the IDs of the local mesh in the context of the
        global parallel mesh. Does not include the IDs of boundary cells.

        E.g., would return [0, 1] for mesh A

        ---------
        | 4 | 5 |     
        ---------  B
        | 2 | 3 |     
        =========
        | 0 | 1 |  A
        ---------
        
        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange((self.offset[1] + self.overlap['bottom']) * self.nx, 
                              (self.offset[1] + self.ny - self.overlap['top']) * self.nx)

    def _getGlobalOverlappingCellIDs(self):
        """
        Return the IDs of the local mesh in the context of the
        global parallel mesh. Includes the IDs of boundary cells.
        
        E.g., would return [0, 1, 2, 3] for mesh A

        ---------
        | 4 | 5 |     
        ---------  B
        | 2 | 3 |     
        =========
        | 0 | 1 |  A
        ---------
        
        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(self.offset[1] * self.nx, (self.offset[1] + self.ny) * self.nx)

    def _getLocalNonOverlappingCellIDs(self):
        """
        Return the IDs of the local mesh in isolation. 
        Does not include the IDs of boundary cells.
        
        E.g., would return [0, 1] for mesh A

        ---------
        | 4 | 5 |     
        ---------  B
        | 2 | 3 |     
        =========
        | 0 | 1 |  A
        ---------
        
        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(self.overlap['bottom'] * self.nx, 
                              (self.ny - self.overlap['top']) * self.nx)

    def _getLocalOverlappingCellIDs(self):
        """
        Return the IDs of the local mesh in isolation. 
        Includes the IDs of boundary cells.
        
        E.g., would return [0, 1, 2, 3] for mesh A

        ---------
        |   |   |     
        ---------  B
        | 2 | 3 |     
        =========
        | 0 | 1 |  A
        ---------
        
        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(0, self.ny * self.nx)
    
## pickling

    def __getstate__(self):
        """
        Used internally to collect the necessary information to ``pickle`` the 
        `Grid2D` to persistent storage.
        """
        return self.args

    def __setstate__(self, dict):
        """
        Used internally to create a new `Grid2D` from ``pickled`` 
        persistent storage.
        """
        self.__init__(**dict)

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
            >>> from fipy.tools import parallel
            >>> print parallel.procID > 0 or numerix.allequal(vertices, mesh._createVertices())
            True
        
            >>> faces = numerix.array(((1, 2, 3, 4, 5, 6, 8, 9, 10, 0, 5, 6, 7, 4, 9, 10, 11),
            ...                        (0, 1, 2, 5, 6, 7, 9, 10, 11, 4, 1, 2, 3, 8, 5, 6, 7)))
            >>> print parallel.procID > 0 or numerix.allequal(faces, mesh._createFaces())
            True

            >>> cells = numerix.array(((0, 1, 2, 3, 4, 5),
            ...                        (10, 11, 12, 14, 15, 16),
            ...                        (3, 4, 5, 6, 7, 8),
            ...                        (9, 10, 11, 13, 14, 15)))
            >>> print parallel.procID > 0 or numerix.allequal(cells, mesh._createCells())
            True

            >>> externalFaces = numerix.array((0, 1, 2, 6, 7, 8, 9 , 12, 13, 16))
            >>> print parallel.procID > 0 or numerix.allequal(externalFaces, 
            ...                                               numerix.nonzero(mesh.getExteriorFaces()))
            True

            >>> internalFaces = numerix.array((3, 4, 5, 10, 11, 14, 15))
            >>> print parallel.procID > 0 or numerix.allequal(internalFaces, 
            ...                                               numerix.nonzero(mesh.getInteriorFaces()))
            True

            >>> from fipy.tools.numerix import MA
            >>> faceCellIds = MA.masked_values(((0, 1, 2, 0, 1, 2, 3, 4, 5, 0, 0, 1, 2, 3, 3, 4, 5),
            ...                                 (-1, -1, -1, 3, 4, 5, -1, -1, -1, -1, 1, 2, -1, -1, 4, 5, -1)), -1)
            >>> print parallel.procID > 0 or numerix.allequal(faceCellIds, mesh.getFaceCellIDs())
            True
            
            >>> faceAreas = numerix.array((dx, dx, dx, dx, dx, dx, dx, dx, dx,
            ...                            dy, dy, dy, dy, dy, dy, dy, dy))
            >>> print parallel.procID > 0 or numerix.allclose(faceAreas, mesh._getFaceAreas(), atol = 1e-10, rtol = 1e-10)
            True
            
            >>> faceCoords = numerix.take(vertices, faces, axis=1)
            >>> faceCenters = (faceCoords[...,0,:] + faceCoords[...,1,:]) / 2.
            >>> print parallel.procID > 0 or numerix.allclose(faceCenters, mesh.getFaceCenters(), atol = 1e-10, rtol = 1e-10)
            True

            >>> faceNormals = numerix.array(((0., 0., 0., 0., 0., 0., 0., 0., 0., -1., 1., 1., 1., -1., 1., 1., 1.),
            ...                              (-1., -1., -1., 1., 1., 1., 1., 1., 1., 0., 0., 0., 0., 0., 0., 0., 0.)))
            >>> print parallel.procID > 0 or numerix.allclose(faceNormals, mesh._getFaceNormals(), atol = 1e-10, rtol = 1e-10)
            True

            >>> cellToFaceOrientations = numerix.array(((1,  1,  1, -1, -1, -1),
            ...                                         (1,  1,  1,  1,  1,  1),
            ...                                         (1,  1,  1,  1,  1,  1),
            ...                                         (1, -1, -1,  1, -1, -1)))
            >>> print parallel.procID > 0 or numerix.allequal(cellToFaceOrientations, mesh._getCellFaceOrientations())
            True
                                             
            >>> cellVolumes = numerix.array((dx*dy, dx*dy, dx*dy, dx*dy, dx*dy, dx*dy))
            >>> print parallel.procID > 0 or numerix.allclose(cellVolumes, mesh.getCellVolumes(), atol = 1e-10, rtol = 1e-10)
            True

            >>> cellCenters = numerix.array(((dx/2., 3.*dx/2., 5.*dx/2., dx/2., 3.*dx/2., 5.*dx/2.),
            ...                              (dy/2., dy/2., dy/2., 3.*dy/2., 3.*dy/2., 3.*dy/2.)))
            >>> print numerix.allclose(mesh.getCellCenters(), cellCenters, atol = 1e-10, rtol = 1e-10)
            True
                                              
            >>> faceToCellDistances = MA.masked_values(((dy / 2., dy / 2., dy / 2., dy / 2., dy / 2., dy / 2., dy / 2., dy / 2., dy / 2., dx / 2., dx / 2., dx / 2., dx / 2., dx / 2., dx / 2., dx / 2., dx / 2.),
            ...                                         (-1, -1, -1, dy / 2., dy / 2., dy / 2., -1, -1, -1, -1, dx / 2., dx / 2., -1, -1, dx / 2., dx / 2., -1)), -1)
            >>> print parallel.procID > 0 or numerix.allclose(faceToCellDistances, mesh._getFaceToCellDistances(), atol = 1e-10, rtol = 1e-10)
            True
                                              
            >>> cellDistances = numerix.array((dy / 2., dy / 2., dy / 2.,
            ...                                dy, dy, dy,
            ...                                dy / 2., dy / 2., dy / 2.,
            ...                                dx / 2., dx, dx,
            ...                                dx / 2.,
            ...                                dx / 2., dx, dx,
            ...                                dx / 2.))
            >>> print parallel.procID > 0 or numerix.allclose(cellDistances, mesh._getCellDistances(), atol = 1e-10, rtol = 1e-10)
            True
            
            >>> faceToCellDistanceRatios = faceToCellDistances[0] / cellDistances
            >>> print parallel.procID > 0 or numerix.allclose(faceToCellDistanceRatios, mesh._getFaceToCellDistanceRatio(), atol = 1e-10, rtol = 1e-10)
            True

            >>> areaProjections = faceNormals * faceAreas
            >>> print parallel.procID > 0 or numerix.allclose(areaProjections, mesh._getAreaProjections(), atol = 1e-10, rtol = 1e-10)
            True

            >>> tangents1 = numerix.array(((1., 1., 1., -1., -1., -1., -1., -1., -1., 0., 0., 0., 0., 0., 0., 0., 0.),
            ...                            (0., 0., 0., 0., 0., 0., 0., 0., 0., -1., 1., 1., 1., -1., 1., 1., 1.)))
            >>> print parallel.procID > 0 or numerix.allclose(tangents1, mesh._getFaceTangents1(), atol = 1e-10, rtol = 1e-10)
            True

            >>> tangents2 = numerix.zeros((2, 17), 'd')
            >>> print parallel.procID > 0 or numerix.allclose(tangents2, mesh._getFaceTangents2(), atol = 1e-10, rtol = 1e-10)
            True

            >>> cellToCellIDs = MA.masked_values(((-1, -1, -1, 0, 1, 2),
            ...                                   (1, 2, -1, 4, 5, -1),
            ...                                   (3, 4, 5, -1, -1, -1),
            ...                                   (-1, 0, 1, -1, 3, 4)), -1)
            >>> print parallel.procID > 0 or numerix.allequal(cellToCellIDs, mesh._getCellToCellIDs())
            True

            >>> cellToCellDistances = MA.masked_values(((dy / 2., dy / 2., dy / 2.,      dy,      dy,      dy),
            ...                                         (     dx,      dx, dx / 2.,      dx,      dx, dx / 2.),
            ...                                         (     dy,      dy,      dy, dy / 2., dy / 2., dy / 2.),
            ...                                         (dx / 2.,      dx,      dx, dx / 2.,      dx,      dx)), -1)
            >>> print parallel.procID > 0 or numerix.allclose(cellToCellDistances, mesh._getCellToCellDistances(), atol = 1e-10, rtol = 1e-10)
            True

            >>> interiorCellIDs = numerix.array(())
            >>> print parallel.procID > 0 or numerix.allequal(interiorCellIDs, mesh._getInteriorCellIDs())
            True

            >>> exteriorCellIDs = numerix.array((0, 1, 2, 3, 4, 5))
            >>> print parallel.procID > 0 or numerix.allequal(exteriorCellIDs, mesh._getExteriorCellIDs())
            True

            >>> cellNormals = numerix.array(((( 0,  0,  0,  0,  0,  0),
            ...                               ( 1,  1,  1,  1,  1,  1),
            ...                               ( 0,  0,  0,  0,  0,  0),
            ...                               (-1, -1, -1, -1, -1, -1)),
            ...                              ((-1, -1, -1, -1, -1, -1),
            ...                               ( 0,  0,  0,  0,  0,  0),
            ...                               ( 1,  1,  1,  1,  1,  1),
            ...                               ( 0,  0,  0,  0,  0,  0))))
            >>> print parallel.procID > 0 or numerix.allclose(cellNormals, mesh._getCellNormals(), atol = 1e-10, rtol = 1e-10)
            True

            >>> cellAreaProjections = numerix.array((((  0,  0,  0,  0,  0,  0),
            ...                                       ( dy, dy, dy, dy, dy, dy),
            ...                                       (  0,  0,  0,  0,  0,  0),
            ...                                       (-dy,-dy,-dy,-dy,-dy,-dy)),
            ...                                      ((-dx,-dx,-dx,-dx,-dx,-dx),
            ...                                       (  0,  0,  0,  0,  0,  0),
            ...                                       ( dx, dx, dx, dx, dx, dx),
            ...                                       (  0,  0,  0,  0,  0,  0))))
            >>> print parallel.procID > 0 or numerix.allclose(cellAreaProjections, mesh._getCellAreaProjections(), atol = 1e-10, rtol = 1e-10)
            True

            >>> cellVertexIDs = MA.masked_values(((5, 6, 7, 9, 10, 11),
            ...                                   (4, 5, 6, 8,  9, 10),
            ...                                   (1, 2, 3, 5,  6,  7),
            ...                                   (0, 1, 2, 4,  5,  6)), -1000)

            >>> print parallel.procID > 0 or numerix.allclose(mesh._getCellVertexIDs(), cellVertexIDs)
            True

            >>> from fipy.tools import dump            
            >>> (f, filename) = dump.write(mesh, extension = '.gz')            
            >>> unpickledMesh = dump.read(filename, f)

            >>> print numerix.allclose(mesh.getCellCenters(), unpickledMesh.getCellCenters())
            True
        """

def _test():
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
