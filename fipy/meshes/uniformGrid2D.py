#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "uniformGrid1D.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #  Author: James O'Beirne <james.obeirne@gmail.com>
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
2D rectangular Mesh with constant spacing in x and constant spacing in y
"""
__docformat__ = 'restructuredtext'

from fipy.meshes.grid2D import Grid2D
from fipy.meshes.topologies import _UniformMeshTopology2D
from fipy.meshes.geometries import _UniformGridGeometry2D

from fipy.tools import numerix
from fipy.tools.numerix import MA
from fipy.tools.dimensions.physicalField import PhysicalField
from fipy.tools import inline
from fipy.tools import parallel
from fipy.tools.decorators import getsetDeprecated

class UniformGrid2D(Grid2D):
    """
    Creates a 2D grid mesh with horizontal faces numbered
    first and then vertical faces.
    """
    def __init__(self, dx=1., dy=1., nx=1, ny=1, origin=((0,),(0,)), 
                 overlap=2, communicator=parallel,
                 GeomClass=_UniformGridGeometry2D):        
        self.args = {
            'dx': dx, 
            'dy': dy, 
            'nx': nx, 
            'ny': ny, 
            'origin': origin,
            'overlap': overlap,
            'communicator': communicator
        }

        self.dim = 2
        
        self.dx = PhysicalField(value = dx)
        scale = PhysicalField(value = 1, unit = self.dx.getUnit())
        self.dx /= scale
        
        nx = int(nx)
        
        self.dy = PhysicalField(value = dy)
        if self.dy.getUnit().isDimensionless():
            self.dy = dy
        else:
            self.dy /= scale
            
        ny = int(ny)

        self.globalNumberOfCells = nx * ny
        self.globalNumberOfFaces = nx * (ny + 1) + ny * (nx + 1)
        
        (self.nx,
         self.ny,
         self.overlap,
         self.offset) = self._calcParallelGridInfo(nx, ny, overlap, communicator)
        
        self.origin = PhysicalField(value = origin)
        self.origin /= scale

        self.origin += ((self.offset[0] * float(self.dx),),
                        (self.offset[1] * float(self.dy),))

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

        self.numberOfHorizontalFaces = self.nx * self.numberOfHorizontalRows
        self.numberOfVerticalFaces = self.numberOfVerticalColumns * self.ny
        self.numberOfFaces = self.numberOfHorizontalFaces + self.numberOfVerticalFaces
        self.numberOfCells = self.nx * self.ny
        
        self._topology = _UniformMeshTopology2D(self)
        self._geometry = GeomClass(self.dx, self.dy,
                                   self.nx, self.ny,
                                   self.origin,
                                   self.numberOfFaces,
                                   self.numberOfHorizontalFaces,
                                   self.numberOfVerticalFaces,
                                   self.numberOfHorizontalRows,
                                   self.numberOfVerticalColumns,
                                   self.numberOfCells)
        
        self.communicator = communicator
        
    def _translate(self, vector):
        return self.__class__(dx = self.args['dx'], nx = self.args['nx'], 
                              dy = self.args['dy'], ny = self.args['ny'], 
                             origin = numerix.array(self.args['origin']) + vector, overlap=self.args['overlap'])

    def __mul__(self, factor):
        if numerix.shape(factor) is ():
            factor = numerix.resize(factor, (2,1))
        
        return UniformGrid2D(dx=self.args['dx'] * numerix.array(factor[0]), nx=self.args['nx'], 
                             dy=self.args['dy'] * numerix.array(factor[1]), ny=self.args['ny'], 
                             origin=numerix.array(self.args['origin']) * factor, overlap=self.args['overlap'])

    def _getConcatenableMesh(self):
        from fipy.meshes.grid2D import Grid2D
        args = self.args.copy()
        origin = args['origin']
        from fipy.tools import serial
        args['communicator'] = serial
        del args['origin']
        return Grid2D(**args) + origin

    @getsetDeprecated
    def _getCellFaceIDs(self):
        return self.cellFaceIDs

    @property
    def cellFaceIDs(self):
        return self._createCells()
        
    @property
    def _maxFacesPerCell(self):
        return 4
        
##         from numMesh/mesh

    def _getVertexCoords(self):
        return self._createVertices() + self.origin

    vertexCoords = property(_getVertexCoords)

    @getsetDeprecated
    def getFaceCellIDs(self):
        return self.faceCellIDs

    @property
    def faceCellIDs(self):
        return inline._optionalInline(self._getFaceCellIDsIn, self._getFaceCellIDsPy)

    def _getFaceCellIDsIn(self):
        faceCellIDs = numerix.zeros((2, self.numberOfFaces))
        mask = numerix.zeros((2, self.numberOfFaces))
        
        inline._runInline("""
            int ID = j * ni + i; 
            int rowlength = ni * nj + Nhor + nj;

            faceCellIDs[ID + 0 * rowlength] = ID - ni;
            faceCellIDs[ID + 1 * rowlength] = ID;

            faceCellIDs[ID + Nhor + j + 0 * rowlength] = ID - 1;
            faceCellIDs[ID + Nhor + j + 1 * rowlength] = ID;

            if (j == 0) {
                faceCellIDs[ID + 0 * rowlength] = ID;
                mask[ID + 1 * rowlength] = 1;
            }

            if (j == nj - 1) {
                faceCellIDs[ID + ni + 0 * rowlength] = ID;
                mask[ID + ni + 1 * rowlength] = 1;
            }

            if (i == 0) {
                faceCellIDs[ID + Nhor + j + 0 * rowlength] = ID;
                mask[ID + Nhor + j + 1 * rowlength] = 1;
            }

            if ( i == ni - 1 ) {
                faceCellIDs[ID + Nhor + j + 1 + 0 * rowlength] = ID;
                mask[ID + Nhor + j + 1 + 1 * rowlength] = 1;
            }
        """,
        Nhor=self.numberOfHorizontalFaces,
        mask=mask,
        faceCellIDs=faceCellIDs,
        ni=self.nx,
        nj=self.ny)

        return MA.masked_where(mask, faceCellIDs)

    def _getFaceCellIDsPy(self):
        Hids = numerix.zeros((2, self.nx, self.numberOfHorizontalRows))
        indices = numerix.indices((self.nx, self.numberOfHorizontalRows))
        Hids[1] = indices[0] + indices[1] * self.nx
        Hids[0] = Hids[1] - self.nx
        if self.numberOfHorizontalRows > 0:
            Hids[0,...,0] = Hids[1,...,0]
            Hids[1,...,0] = -1
            Hids[1,...,-1] = -1

        Vids = numerix.zeros((2, self.numberOfVerticalColumns, self.ny))
        indices = numerix.indices((self.numberOfVerticalColumns, self.ny))
        Vids[1] = indices[0] + indices[1] * self.nx
        Vids[0] = Vids[1] - 1
        if self.numberOfVerticalColumns > 0:
            Vids[0,0] = Vids[1,0]
            Vids[1,0] = -1
            Vids[1,-1] = -1
        
        return MA.masked_values(numerix.concatenate((Hids.reshape((2, self.numberOfHorizontalFaces), order="FORTRAN"), 
                                                     Vids.reshape((2, self.numberOfFaces - self.numberOfHorizontalFaces), order="FORTRAN")), axis=1), value = -1)
    
    @getsetDeprecated
    def _getCellVertexIDs(self):
        return self._cellVertexIDs

    @property
    def _cellVertexIDs(self):
        ids = numerix.zeros((4, self.nx, self.ny))
        indices = numerix.indices((self.nx, self.ny))
        ids[1] = indices[0] + (indices[1] + 1) * self.numberOfVerticalColumns
        ids[0] = ids[1] + 1
        ids[3] = indices[0] + indices[1] * self.numberOfVerticalColumns
        ids[2] = ids[3] + 1
        
        return numerix.reshape(ids, (4, self.numberOfCells))
        
    @getsetDeprecated
    def _getFaceVertexIDs(self):
        return self.faceVertexIDs

    @property
    def faceVertexIDs(self):
        Hids = numerix.zeros((2, self.nx, self.numberOfHorizontalRows))
        indices = numerix.indices((self.nx, self.numberOfHorizontalRows))
        Hids[0] = indices[0] + indices[1] * self.numberOfVerticalColumns
        Hids[1] = Hids[0] + 1

        Vids = numerix.zeros((2, self.numberOfVerticalColumns, self.ny))
        indices = numerix.indices((self.numberOfVerticalColumns, self.ny))
        Vids[0] = indices[0] + indices[1] * self.numberOfVerticalColumns
        Vids[1] = Vids[0] + self.numberOfVerticalColumns
        
        return numerix.concatenate((Hids.reshape((2, self.numberOfHorizontalFaces), order="FORTRAN"), 
                                    Vids.reshape((2, self.numberOfFaces - self.numberOfHorizontalFaces),
                                                 order="FORTRAN")),
                                   axis=1)

    @getsetDeprecated
    def _getOrderedCellVertexIDs(self):
        return self._orderedCellVertexIDs

    @property
    def _orderedCellVertexIDs(self):
        ids = numerix.zeros((4, self.nx, self.ny))
        indices = numerix.indices((self.nx, self.ny))
        ids[2] = indices[0] + (indices[1] + 1) * self.numberOfVerticalColumns
        ids[1] = ids[2] + 1
        ids[3] = indices[0] + indices[1] * self.numberOfVerticalColumns
        ids[0] = ids[3] + 1
        
        return ids.reshape((4, self.numberOfCells), order="FORTRAN")
        
##     scaling
    
    def _setScaledGeometry(self):
        pass

    def _getNearestCellID(self, points):
        """
        Test cases

           >>> from fipy import *
           >>> m = Grid2D(nx=3, ny=2)
           >>> eps = numerix.array([[1e-5, 1e-5]])
           >>> print m._getNearestCellID(((0., .9, 3.), (0., 2., 2.)))
           [0 3 5]
           >>> print m._getNearestCellID(([1.1], [1.5]))
           [4]
           >>> m0 = Grid2D(nx=2, ny=2, dx=1., dy=1.)
           >>> m1 = Grid2D(nx=4, ny=4, dx=.5, dy=.5)
           >>> print m0._getNearestCellID(m1.getCellCenters().getGlobalValue())
           [0 0 1 1 0 0 1 1 2 2 3 3 2 2 3 3]
           
        """
        nx = self.args['nx']
        ny = self.args['ny']
        
        if nx == 0 or ny == 0:
            return numerix.arange(0)
            
        x0, y0 = self.getCellCenters().getGlobalValue()[...,0]        
        xi, yi = points
        dx, dy = self.dx, self.dy
        
        i = numerix.array(numerix.rint(((xi - x0) / dx)), 'l')
        i[i < 0] = 0
        i[i > nx - 1] = nx - 1

        j = numerix.array(numerix.rint(((yi - y0) / dy)), 'l')
        j[j < 0] = 0
        j[j > ny - 1]  = ny - 1

        return j * nx + i

    def _test(self):
        """
        These tests are not useful as documentation, but are here to ensure
        everything works as expected.

            >>> dx = 0.5
            >>> dy = 2.
            >>> nx = 3
            >>> ny = 2
            
            >>> mesh = UniformGrid2D(nx = nx, ny = ny, dx = dx, dy = dy)     
            
            >>> vertices = numerix.array(((0., 1., 2., 3., 0., 1., 2., 3., 0., 1., 2., 3.),
            ...                           (0., 0., 0., 0., 1., 1., 1., 1., 2., 2., 2., 2.)))
            >>> vertices *= numerix.array(((dx,), (dy,)))
            >>> from fipy.tools import parallel
            >>> print parallel.procID > 0 or numerix.allequal(vertices, mesh._createVertices())
            True
        
            >>> faces = numerix.array(((1, 2, 3, 4, 5, 6, 8, 9, 10, 0, 5, 6, 7, 4, 9, 10, 11),
            ...                        (0, 1, 2, 5, 6, 7, 9, 10, 11, 4, 1, 2, 3, 8, 5, 6, 7)))
            >>> print parallel.procID > 0 or numerix.allequal(faces, mesh._createFaces()[0])
            True

            >>> cells = numerix.array(((0, 1, 2, 3, 4, 5),
            ...                        (10, 11, 12, 14, 15, 16),
            ...                        (3, 4, 5, 6, 7, 8),
            ...                        (9, 10, 11, 13, 14, 15)))
            >>> print parallel.procID > 0 or numerix.allequal(cells, mesh._createCells())
            True

            >>> externalFaces = numerix.array((0, 1, 2, 6, 7, 8, 9 , 12, 13, 16))
            >>> print parallel.procID > 0 or numerix.allequal(externalFaces, 
            ...                                               numerix.nonzero(mesh.exteriorFaces))
            True

            >>> internalFaces = numerix.array((3, 4, 5, 10, 11, 14, 15))
            >>> print parallel.procID > 0 or numerix.allequal(internalFaces, 
            ...                                               numerix.nonzero(mesh.interiorFaces))
            True

            >>> from fipy.tools.numerix import MA
            >>> faceCellIds = MA.masked_values(((0, 1, 2, 0, 1, 2, 3, 4, 5, 0, 0, 1, 2, 3, 3, 4, 5),
            ...                                 (-1, -1, -1, 3, 4, 5, -1, -1, -1, -1, 1, 2, -1, -1, 4, 5, -1)), -1)
            >>> print parallel.procID > 0 or numerix.allequal(faceCellIds, mesh.faceCellIDs)
            True
            
            >>> faceAreas = numerix.array((dx, dx, dx, dx, dx, dx, dx, dx, dx,
            ...                            dy, dy, dy, dy, dy, dy, dy, dy))
            >>> print parallel.procID > 0 or numerix.allclose(faceAreas, mesh._getFaceAreas(), atol = 1e-10, rtol = 1e-10)
            True
            
            >>> faceCoords = numerix.take(vertices, faces, axis=1)
            >>> faceCenters = (faceCoords[...,0,:] + faceCoords[...,1,:]) / 2.
            >>> print parallel.procID > 0 or numerix.allclose(faceCenters, mesh.faceCenters, atol = 1e-10, rtol = 1e-10)
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
            >>> print parallel.procID > 0 or numerix.allclose(cellVolumes, mesh.cellVolumes, atol = 1e-10, rtol = 1e-10)
            True

            >>> cellCenters = numerix.array(((dx/2., 3.*dx/2., 5.*dx/2., dx/2., 3.*dx/2., 5.*dx/2.),
            ...                              (dy/2., dy/2., dy/2., 3.*dy/2., 3.*dy/2., 3.*dy/2.)))
            >>> print numerix.allclose(cellCenters, mesh.cellCenters, atol = 1e-10, rtol = 1e-10)
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

            >>> cellVertexIDs = MA.masked_array(((5, 6, 7, 9, 10, 11),
            ...                                  (4, 5, 6, 8, 9, 10),
            ...                                  (1, 2, 3, 5, 6, 7),
            ...                                  (0, 1, 2, 4, 5, 6)), -1000)

            >>> print parallel.procID > 0 or numerix.allclose(mesh._getCellVertexIDs(), cellVertexIDs)
            True

            >>> from fipy.tools import dump            
            >>> (f, filename) = dump.write(mesh, extension = '.gz')            
            >>> unpickledMesh = dump.read(filename, f)

            >>> print numerix.allequal(mesh.cellCenters, unpickledMesh.getCellCenters())
            True
            
            >>> faceVertexIDs = [[ 0, 1, 2, 4, 5, 6, 8, 9, 10, 0, 1, 2, 3, 4, 5, 6, 7],
            ...                  [ 1, 2, 3, 5, 6, 7, 9, 10, 11, 4, 5, 6, 7, 8, 9, 10, 11]]
            >>> print parallel.procID > 0 or numerix.allequal(mesh.faceVertexIDs, faceVertexIDs)
            True

            >>> mesh = UniformGrid2D(nx=3)
            >>> print parallel.procID > 0 or numerix.allequal(mesh._getAdjacentCellIDs()[0],
            ...                                               [0, 1, 2, 0, 1, 2, 0, 0, 1, 2])
            True
            >>> print parallel.procID > 0 or numerix.allequal(mesh._getAdjacentCellIDs()[1],
            ...                                               [0, 1, 2, 0, 1, 2, 0, 1, 2, 2])
            True
            
            >>> faceCellIDs = [[0, 1, 2, 0, 1, 2, 0, 0, 1, 2],
            ...                [-1, -1, -1, -1, -1, -1, -1, 1, 2, -1]]
            >>> print parallel.procID > 0 or numerix.allequal(mesh.faceCellIDs.filled(-1),
            ...                                               faceCellIDs)
            True

            >>> mesh = UniformGrid2D(ny=3)
            >>> print parallel.procID > 0 or numerix.allequal(mesh._getAdjacentCellIDs()[0],
            ...                                               [0, 0, 1, 2, 0, 0, 1, 1, 2, 2])
            True
            >>> print parallel.procID > 0 or numerix.allequal(mesh._getAdjacentCellIDs()[1],
            ...                                               [0, 1, 2, 2, 0, 0, 1, 1, 2, 2])
            True
            >>> faceCellIDs = [[0, 0, 1, 2, 0, 0, 1, 1, 2, 2],
            ...                [-1, 1, 2, -1, -1, -1, -1, -1, -1, -1]]
            >>> print parallel.procID > 0 or numerix.allequal(mesh.faceCellIDs.filled(-1),
            ...                                               faceCellIDs)
            True

        Following test added to change nx, ny argment to integer when its a float to prevent
        warnings from the solver.

            >>> from fipy import *
            >>> mesh = UniformGrid2D(nx=3., ny=3., dx=1., dy=1.)
            >>> var = CellVariable(mesh=mesh)
            >>> DiffusionTerm().solve(var)

        """

def _test():
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
