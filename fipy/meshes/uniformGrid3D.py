#!/usr/bin/env python

## -*-Pyth-*-
 # ########################################################################
 # FiPy - a finite volume PDE solver in Python
 # 
 # FILE: "uniformGrid3D.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #  Author: James O'Beirne <james.obeirne@gmail.com>
 #   mail: NIST
 #    www: <http://www.ctcms.nist.gov/fipy/>
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
 # ########################################################################
 ##

from fipy.tools.numerix import MA

from fipy.meshes.grid3D import Grid3D
from fipy.meshes.topologies import _UniformMeshTopology3D
from fipy.meshes.geometries import _UniformGridGeometry3D
from fipy.tools import numerix
from fipy.tools.dimensions.physicalField import PhysicalField
from fipy.tools.decorators import getsetDeprecated

from fipy.tools import parallel

class UniformGrid3D(Grid3D):
    """
    3D rectangular-prism Mesh with uniform grid spacing in each dimension.

    X axis runs from left to right.
    Y axis runs from bottom to top.
    Z axis runs from front to back.

    Numbering System:

    Vertices: Numbered in the usual way. X coordinate changes most quickly, then Y, then Z.
    
    *** arrays are arranged Z, Y, X because in numerix, the final index is the one that changes the most quickly ***

    Cells: Same numbering system as vertices.

    Faces: XY faces numbered first, then XZ faces, then YZ faces. Within each subcategory, it is numbered in the usual way.
    """
    def __init__(self, dx = 1., dy = 1., dz = 1., nx = 1, ny = 1, nz = 1, 
                 origin = [[0], [0], [0]], overlap=2, communicator=parallel):
        self.args = {
            'dx': dx, 
            'dy': dy,
            'dz': dz,
            'nx': nx, 
            'ny': ny,
            'nz': nz,
            'origin': origin,
            'overlap': overlap,
            'communicator': communicator
        }
        
        self.dim = 3
        
        self.dx = PhysicalField(value = dx)
        scale = PhysicalField(value = 1, unit = self.dx.unit)
        self.dx /= scale
        
        nx = int(nx)
        
        self.dy = PhysicalField(value = dy)
        if self.dy.unit.isDimensionless():
            self.dy = dy
        else:
            self.dy /= scale
            
        ny = int(ny)
        
        self.dz = PhysicalField(value = dy)
        if self.dz.unit.isDimensionless():
            self.dz = dz
        else:
            self.dz /= scale
            
        nz = int(nz)

        self.globalNumberOfCells = nx * ny * nz
        self.globalNumberOfFaces = nx * nz * (ny + 1) + ny * nz * (nx + 1) \
                                     + nx * ny * (nz + 1)

        (self.nx,
         self.ny,
         self.nz,
         self.overlap,
         self.offset) = self._calcParallelGridInfo(nx, ny, nz, overlap, communicator)
        
        self.origin = PhysicalField(value = origin)
        self.origin /= scale

        self.origin += ((self.offset[0] * float(self.dx),),
                        (self.offset[1] * float(self.dy),),
                        (self.offset[2] * float(self.dz),))

        if self.nx == 0 or self.ny == 0 or self.nz == 0:
            self.nx = 0
            self.ny = 0
            self.nz = 0
        if self.nx == 0 or self.ny == 0 or self.nz == 0:
            self.numberOfHorizontalRows = 0
            self.numberOfVerticalColumns = 0
            self.numberOflayers = 0
        else:
            self.numberOfHorizontalRows = (self.ny + 1)
            self.numberOfVerticalColumns = (self.nx + 1)
            self.numberOfLayers = (self.nz + 1)

        self.numberOfVertices = (self.nx + 1) * (self.ny + 1) * (self.nz + 1)
        self.numberOfXYFaces = self.nx * self.ny * (self.nz + 1)
        self.numberOfXZFaces = self.nx * (self.ny + 1) * self.nz
        self.numberOfYZFaces = (self.nx + 1) * self.ny * self.nz
        self.numberOfFaces = self.numberOfXYFaces + self.numberOfXZFaces + self.numberOfYZFaces
        self.numberOfCells = self.nx * self.ny * self.nz

        self._topology = _UniformMeshTopology3D(self.nx, self.ny, self.nz,
                                                self.numberOfCells,
                                                self._maxFacesPerCell,
                                                self._XYFaceIDs,
                                                self._XZFaceIDs,
                                                self._YZFaceIDs,
                                                self.faceCellIDs,
                                                self.cellFaceIDs,
                                                self)

        self._geometry = _UniformGridGeometry3D(self.dx, self.dy, self.dz,
                                               self.nx, self.ny, self.nz,
                                               self.numberOfCells,
                                               self.numberOfXYFaces,
                                               self.numberOfXZFaces,
                                               self.numberOfYZFaces,
                                               self.origin)
        
        self.communicator = communicator
        
    def _translate(self, vector):
        return self.__class__(dx = self.args['dx'], nx = self.args['nx'], 
                              dy = self.args['dy'], ny = self.args['ny'],
                              dz = self.args['dz'], nz = self.args['nz'],
                             origin = numerix.array(self.args['origin']) + vector, overlap=self.args['overlap'])
        
    def __mul__(self, factor):
        return UniformGrid3D(dx = self.dx * factor, nx = self.nx, 
                             dy = self.dy * factor, ny = self.ny, 
                             dz = self.dz * factor, nz = self.nz, 
                             origin = self.origin * factor)

    @property
    def _concatenableMesh(self):
        from fipy.meshes.grid3D import Grid3D
        args = self.args.copy()
        origin = args['origin']
        from fipy.tools import serial
        args['communicator'] = serial
        del args['origin']
        return Grid3D(**args) + origin

##     get topology methods

##         from common/mesh
        
    @getsetDeprecated(new_name="cellFaceIDs")
    def _getCellFaceIDs(self):
        return self.cellFaceIDs

    @property
    def cellFaceIDs(self):
        return MA.array(self._createCells())

    @getsetDeprecated
    def _getXYFaceIDs(self):
        return self._XYFaceIDs

    @property
    def _XYFaceIDs(self):
        ids = numerix.arange(0, self.numberOfXYFaces)
        return ids.reshape((self.nz + 1, self.ny, self.nx)).swapaxes(0,2)

    @getsetDeprecated
    def _getXZFaceIDs(self):
        return self._XZFaceIDs

    @property
    def _XZFaceIDs(self):
        ids = numerix.arange(self.numberOfXYFaces, self.numberOfXYFaces + self.numberOfXZFaces)
        return ids.reshape((self.nz, self.ny + 1, self.nx)).swapaxes(0,2)

    @getsetDeprecated
    def _getYZFaceIDs(self):
        return self._YZFaceIDs

    @property
    def _YZFaceIDs(self):
        ids = numerix.arange(self.numberOfXYFaces + self.numberOfXZFaces, self.numberOfFaces)
        return ids.reshape((self.nz, self.ny, self.nx + 1)).swapaxes(0,2)
   
    @property
    def _maxFacesPerCell(self):
        return 6
        
##         from numMesh/mesh

    @getsetDeprecated
    def _getVertexCoords(self):
        return self.vertexCoords

    @property
    def vertexCoords(self):
        return self._createVertices() + self.origin

    @getsetDeprecated
    def getFaceCellIDs(self):
        return self.faceCellIDs

    @property
    def faceCellIDs(self):
        XYids = MA.zeros((2, self.nx, self.ny, self.nz + 1), 'l')
        indices = numerix.indices((self.nx, self.ny, self.nz + 1))
        XYids[1] = indices[0] + (indices[1] + indices[2] * self.ny) * self.nx
        XYids[0] = XYids[1] - self.nx * self.ny
        XYids[0,..., 0] = XYids[1,..., 0]
        XYids[1,..., 0] = MA.masked
        XYids[1,...,-1] = MA.masked
        
        XZids = MA.zeros((2, self.nx, self.ny + 1, self.nz), 'l')
        indices = numerix.indices((self.nx, self.ny + 1, self.nz))
        XZids[1] = indices[0] + (indices[1] + indices[2] * self.ny) * self.nx
        XZids[0] = XZids[1] - self.nx
        XZids[0,..., 0,...] = XZids[1,..., 0,...]
        XZids[1,..., 0,...] = MA.masked
        XZids[1,...,-1,...] = MA.masked

        YZids = MA.zeros((2, self.nx + 1, self.ny, self.nz), 'l')
        indices = numerix.indices((self.nx + 1, self.ny, self.nz))
        YZids[1] = indices[0] + (indices[1] + indices[2] * self.ny) * self.nx
        YZids[0] = YZids[1] - 1
        YZids[0, 0] = YZids[1, 0]
        YZids[1, 0] = MA.masked
        YZids[1,-1] = MA.masked

        return MA.concatenate((XYids.swapaxes(1,3).reshape((2, self.numberOfXYFaces)), 
                               XZids.swapaxes(1,3).reshape((2, self.numberOfXZFaces)), 
                               YZids.swapaxes(1,3).reshape((2, self.numberOfYZFaces))), axis=1)

##         from common/mesh
                                   
    @getsetDeprecated
    def _getCellVertexIDs(self):
        return self._cellVertexIDs

    @property
    def _cellVertexIDs(self):
        ids = numerix.zeros((8, self.nx, self.ny, self.nz))
        indices = numerix.indices((self.nx, self.ny, self.nz))
        ids[1] = indices[0] + (indices[1] + (indices[2] + 1) * (self.ny + 1) + 1) * (self.nx + 1)
        ids[0] = ids[1] + 1
        ids[3] = indices[0] + (indices[1] + (indices[2] + 1) * (self.ny + 1)) * (self.nx + 1)
        ids[2] = ids[3] + 1
        ids[5] = indices[0] + (indices[1] + indices[2] * (self.ny + 1) + 1) * (self.nx + 1)
        ids[4] = ids[5] + 1
        ids[7] = indices[0] + (indices[1] + indices[2] * (self.ny + 1)) * (self.nx + 1)
        ids[6] = ids[7] + 1
        
        return numerix.reshape(ids.swapaxes(1,3), (8, self.numberOfCells))
        
    @getsetDeprecated(new_name="faceVertexIDs")
    def _getFaceVertexIDs(self):
        return self.faceVertexIDs

    @property
    def faceVertexIDs(self):
       return self._createFaces()[1]

    @property
    def _orderedCellVertexIDs(self):
        """Correct ordering for VTK_VOXEL"""
        return self._cellVertexIDs     
        
##     scaling
    
    def _setScaledGeometry(self):
        pass
    
    def _getNearestCellID(self, points):
        nx = self.args['nx']
        ny = self.args['ny']
        nz = self.args['nz']
        
        x0, y0, z0 = self.cellCenters[...,0]        
        xi, yi, zi = points
        nx, ny, nz = self.shape
        dx, dy, dz = self.dx, self.dy, self.dz
        
        i = numerix.array(numerix.rint(((xi - x0) / dx)), 'l')
        i[i < 0] = 0
        i[i > nx - 1] = nx - 1

        j = numerix.array(numerix.rint(((yi - y0) / dy)), 'l')
        j[j < 0] = 0
        j[j > ny - 1]  = ny - 1

        k = numerix.array(numerix.rint(((zi - z0) / dz)), 'l')
        k[k < 0] = 0
        k[k > nz - 1]  = nz - 1
        
        return k * ny * nx + j * nx + i
        
    def _test(self):
        """
        These tests are not useful as documentation, but are here to ensure
        everything works as expected.
        
            >>> dx = 0.5
            >>> dy = 2.
            >>> dz = 4.
            >>> nx = 3
            >>> ny = 2
            >>> nz = 1

            >>> mesh = UniformGrid3D(nx=nx, ny=ny, nz=nz, dx=dx, dy=dy, dz=dz)

            >>> XYFaceIDs = numerix.array((((0, 6), (3, 9)),
            ...                            ((1, 7), (4, 10)),
            ...                            ((2, 8), (5, 11))))
            >>> print parallel.procID > 0 or numerix.allequal(XYFaceIDs, mesh._XYFaceIDs)
            True
              
            >>> XZFaceIDs = numerix.array((((12,), (15,), (18,)),
            ...                            ((13,), (16,), (19,)),
            ...                            ((14,), (17,), (20,))))
            >>> print parallel.procID > 0 or numerix.allequal(mesh._XZFaceIDs, XZFaceIDs)
            True
            
            >>> YZFaceIDs = numerix.array((((21,), (25,)),
            ...                            ((22,), (26,)),
            ...                            ((23,), (27,)),
            ...                            ((24,), (28,))))
            >>> print parallel.procID > 0 or numerix.allequal(mesh._YZFaceIDs, YZFaceIDs)
            True

            >>> adjacentCellIDs = (numerix.array([0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 0, 1, 2, 3, 4, 5, 0, 0,
            ... 1, 2, 3, 3, 4, 5]), numerix.array([0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 3, 4, 5, 0, 1,
            ... 2, 2, 3, 4, 5, 5]))
            >>> print parallel.procID > 0 or numerix.allequal(mesh._adjacentCellIDs, adjacentCellIDs)
            True

            >>> vertices = numerix.array(((0., 1., 2., 3., 0., 1., 2., 3., 0., 1., 2., 3., 0., 1., 2., 3., 0., 1., 2., 3., 0., 1., 2., 3.),
            ...                           (0., 0., 0., 0., 1., 1., 1., 1., 2., 2., 2., 2., 0., 0., 0., 0., 1., 1., 1., 1., 2., 2., 2., 2.),
            ...                           (0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.)))
            >>> vertices *= numerix.array([[dx], [dy], [dz]])
            
            >>> print parallel.procID > 0 or numerix.allequal(vertices, mesh._createVertices())
            True
        
            >>> faces = numerix.array((( 0,  1,  2,  4,  5,  6, 12, 13, 14, 16, 17, 18,  0,  1,  2,  4,  5,  6,  8,  9, 10,  0,  1,  2,  3,  4,  5,  6,  7),
            ...                        ( 1,  2,  3,  5,  6,  7, 13, 14, 15, 17, 18, 19,  1,  2,  3,  5,  6,  7,  9, 10, 11,  4,  5,  6,  7,  8,  9, 10, 11),
            ...                        ( 5,  6,  7,  9, 10, 11, 17, 18, 19, 21, 22, 23, 13, 14, 15, 17, 18, 19, 21, 22, 23, 16, 17, 18, 19, 20, 21, 22, 23),
            ...                        ( 4,  5,  6,  8,  9, 10, 16, 17, 18, 20, 21, 22, 12, 13, 14, 16, 17, 18, 20, 21, 22, 12, 13, 14, 15, 16, 17, 18, 19))) 
            >>> print parallel.procID > 0 or numerix.allequal(faces, mesh._createFaces()[1])
            True

            >>> cells = numerix.array(((21, 22, 23, 25, 26, 27),
            ...                        (22, 23, 24, 26, 27, 28),
            ...                        (12, 13, 14, 15, 16, 17),
            ...                        (15, 16, 17, 18, 19, 20),
            ...                        ( 0,  1,  2,  3,  4,  5),
            ...                        ( 6,  7,  8,  9, 10, 11)))
            >>> print parallel.procID > 0 or numerix.allequal(cells, mesh._createCells())
            True

            >>> externalFaces = numerix.array((0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 18, 19, 20, 21, 24, 25, 28))
            >>> print parallel.procID > 0 or numerix.allequal(externalFaces, 
            ...                              numerix.nonzero(mesh.exteriorFaces))
            True

            >>> internalFaces = numerix.array((15, 16, 17, 22, 23, 26, 27))
            >>> print parallel.procID > 0 or numerix.allequal(internalFaces, 
            ...                              numerix.nonzero(mesh.interiorFaces))
            True

            >>> from fipy.tools.numerix import MA
            >>> faceCellIds = MA.masked_values((( 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 0, 1, 2, 3, 4, 5, 0, 0, 1, 2, 3, 3, 4, 5),
            ...                                 (-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 3, 4, 5,-1,-1,-1,-1, 1, 2,-1,-1, 4, 5,-1)), -1) 
            >>> print parallel.procID > 0 or numerix.allequal(faceCellIds, mesh.faceCellIDs)
            True
            
            >>> xy = dx * dy
            >>> xz = dx * dz
            >>> yz = dy * dz
            >>> faceAreas = numerix.array((xy, xy, xy, xy, xy, xy, xy, xy, xy, xy, xy, xy,
            ...                            xz, xz, xz, xz, xz, xz, xz, xz, xz,
            ...                            yz, yz, yz, yz, yz, yz, yz, yz))
            >>> print parallel.procID > 0 or numerix.allclose(faceAreas, mesh._faceAreas, atol = 1e-10, rtol = 1e-10)
            1
            
            >>> faceCoords = numerix.take(vertices, faces, axis=1)
            >>> faceCenters = (faceCoords[...,0,:] + faceCoords[...,1,:] + faceCoords[...,2,:] + faceCoords[...,3,:]) / 4.
            >>> print parallel.procID > 0 or numerix.allclose(faceCenters, mesh.faceCenters, atol = 1e-10, rtol = 1e-10)
            True

            >>> faceNormals = numerix.array((( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 1, 1,-1, 1, 1, 1),
            ...                              ( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0),
            ...                              (-1,-1,-1,-1,-1,-1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)))
            >>> print parallel.procID > 0 or numerix.allclose(faceNormals, mesh._faceNormals, atol = 1e-10, rtol = 1e-10)
            True

            >>> cellToFaceOrientations = numerix.array((( 1,-1,-1, 1,-1,-1),
            ...                                         ( 1, 1, 1, 1, 1, 1),
            ...                                         ( 1, 1, 1,-1,-1,-1),
            ...                                         ( 1, 1, 1, 1, 1, 1),
            ...                                         ( 1, 1, 1, 1, 1, 1),
            ...                                         ( 1, 1, 1, 1, 1, 1)))
            >>> print parallel.procID > 0 or numerix.allequal(cellToFaceOrientations, mesh._cellToFaceOrientations)
            True
                                             
            >>> cellVolumes = numerix.array((dx*dy*dz, dx*dy*dz, dx*dy*dz, dx*dy*dz, dx*dy*dz, dx*dy*dz))
            >>> print parallel.procID > 0 or numerix.allclose(cellVolumes, mesh.cellVolumes, atol = 1e-10, rtol = 1e-10)
            True

            >>> cellCenters = numerix.array(((   dx/2., 3.*dx/2., 5.*dx/2.,   dx/2., 3.*dx/2., 5.*dx/2.),
            ...                              (   dy/2.,    dy/2.,    dy/2.,3.*dy/2., 3.*dy/2., 3.*dy/2.),
            ...                              (   dz/2.,    dz/2.,    dz/2.,   dz/2.,    dz/2.,    dz/2.)))
            >>> print numerix.allclose(cellCenters, mesh.cellCenters, atol = 1e-10, rtol = 1e-10)
            True
                                              
            >>> cellDistances = numerix.array((dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2,
            ...                                dy/2, dy/2, dy/2, dy, dy, dy, dy/2, dy/2, dy/2,
            ...                                dx/2, dx, dx, dx/2, dx/2, dx, dx, dx/2))
            >>> print parallel.procID > 0 or numerix.allclose(cellDistances, mesh._cellDistances, atol = 1e-10, rtol = 1e-10)
            True
            
            >>> faceToCellDistances = MA.masked_values(((dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dy/2, dy/2, dy/2, dy/2, dy/2, dy/2, dy/2, dy/2, dy/2, dx/2, dx/2, dx/2, dx/2, dx/2, dx/2, dx/2, dx/2),
            ...                                         (  -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1, dy/2, dy/2, dy/2,   -1,   -1,   -1,   -1, dx/2, dx/2,   -1,   -1, dx/2, dx/2,   -1)), -1)
            >>> faceToCellDistanceRatios = faceToCellDistances[0] / cellDistances
            >>> print parallel.procID > 0 or numerix.allclose(faceToCellDistanceRatios, mesh._faceToCellDistanceRatio, atol = 1e-10, rtol = 1e-10)
            True

            >>> areaProjections = faceNormals * faceAreas
            >>> print parallel.procID > 0 or numerix.allclose(areaProjections, mesh._areaProjections, atol = 1e-10, rtol = 1e-10)
            True

            >>> tangents1 = numerix.array(((1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0),
            ...                            (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1),
            ...                            (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)))
            >>> print parallel.procID > 0 or numerix.allclose(tangents1, mesh._faceTangents1, atol = 1e-10, rtol = 1e-10)
            True

            >>> tangents2 = numerix.array(((0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
            ...                            (1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
            ...                            (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)))
            >>> print parallel.procID > 0 or numerix.allclose(tangents2, mesh._faceTangents2, atol = 1e-10, rtol = 1e-10)
            True

            >>> cellToCellIDs = MA.masked_values(((-1, 0, 1, -1, 3, 4),
            ...                                   (1, 2, -1, 4, 5, -1),
            ...                                   (-1, -1, -1, 0, 1, 2),
            ...                                   (3, 4, 5, -1, -1, -1),
            ...                                   (-1, -1, -1, -1, -1, -1),
            ...                                   (-1, -1, -1, -1, -1, -1)), -1)
            >>> print parallel.procID > 0 or numerix.allequal(mesh._cellToCellIDs, cellToCellIDs)
            True

            >>> cellToCellIDsFilled = numerix.array([[0, 0, 1, 3, 3, 4],
            ...                                      [1, 2, 2, 4, 5, 5],
            ...                                      [0, 1, 2, 0, 1, 2],
            ...                                      [3, 4, 5, 3, 4, 5],
            ...                                      [0, 1, 2, 3, 4, 5],
            ...                                      [0, 1, 2, 3, 4, 5]])
            >>> print parallel.procID > 0 or numerix.allequal(mesh._cellToCellIDsFilled, cellToCellIDsFilled)
            True
              
            >>> cellToCellDistances = numerix.take(cellDistances, cells)
            >>> print parallel.procID > 0 or numerix.allclose(cellToCellDistances, mesh._cellToCellDistances, atol = 1e-10, rtol = 1e-10)
            True

            >>> cellNormals = numerix.array((((-1, -1, -1, -1, -1, -1),
            ...                               ( 1,  1,  1,  1,  1,  1),
            ...                               ( 0,  0,  0,  0,  0,  0),
            ...                               ( 0,  0,  0,  0,  0,  0),
            ...                               ( 0,  0,  0,  0,  0,  0),
            ...                               ( 0,  0,  0,  0,  0,  0)),
            ...                              (( 0,  0,  0,  0,  0,  0),
            ...                               ( 0,  0,  0,  0,  0,  0),
            ...                               (-1, -1, -1, -1, -1, -1),
            ...                               ( 1,  1,  1,  1,  1,  1),
            ...                               ( 0,  0,  0,  0,  0,  0),
            ...                               ( 0,  0,  0,  0,  0,  0)),
            ...                              (( 0,  0,  0,  0,  0,  0),
            ...                               ( 0,  0,  0,  0,  0,  0),
            ...                               ( 0,  0,  0,  0,  0,  0),
            ...                               ( 0,  0,  0,  0,  0,  0),
            ...                               (-1, -1, -1, -1, -1, -1),
            ...                               ( 1,  1,  1,  1,  1,  1))))
            >>> print parallel.procID > 0 or numerix.allclose(cellNormals, mesh._cellNormals, atol = 1e-10, rtol = 1e-10)
            True

            >>> cellAreaProjections = numerix.array((((-yz,-yz,-yz,-yz,-yz,-yz),
            ...                                       ( yz, yz, yz, yz, yz, yz),
            ...                                       (  0,  0,  0,  0,  0,  0),
            ...                                       (  0,  0,  0,  0,  0,  0),
            ...                                       (  0,  0,  0,  0,  0,  0),
            ...                                       (  0,  0,  0,  0,  0,  0)),
            ...                                      ((  0,  0,  0,  0,  0,  0),
            ...                                       (  0,  0,  0,  0,  0,  0),
            ...                                       (-xz,-xz,-xz,-xz,-xz,-xz),
            ...                                       ( xz, xz, xz, xz, xz, xz),
            ...                                       (  0,  0,  0,  0,  0,  0),
            ...                                       (  0,  0,  0,  0,  0,  0)),
            ...                                      ((  0,  0,  0,  0,  0,  0),
            ...                                       (  0,  0,  0,  0,  0,  0),
            ...                                       (  0,  0,  0,  0,  0,  0),
            ...                                       (  0,  0,  0,  0,  0,  0),
            ...                                       (-xy,-xy,-xy,-xy,-xy,-xy),
            ...                                       ( xy, xy, xy, xy, xy, xy))))
            >>> print parallel.procID > 0 or numerix.allclose(cellAreaProjections, mesh._cellAreaProjections, atol = 1e-10, rtol = 1e-10)
            True

            >>> cellVertexIDs = numerix.array((17, 16, 13, 12, 5, 4, 1, 0))
            >>> cellVertexIDs = numerix.array((cellVertexIDs, cellVertexIDs + 1, cellVertexIDs + 2,
            ...                                cellVertexIDs + 4, cellVertexIDs + 5, cellVertexIDs + 6))
            >>> cellVertexIDs = cellVertexIDs.swapaxes(0,1)
            >>> print parallel.procID > 0 or numerix.allclose(mesh._cellVertexIDs, cellVertexIDs)
            True

            >>> from fipy.tools import dump            
            >>> (f, filename) = dump.write(mesh, extension = '.gz')            
            >>> unpickledMesh = dump.read(filename, f)

            >>> print numerix.allequal(mesh.cellCenters, unpickledMesh.cellCenters)
            True
            
            # Bug #130 & #135 are because we only checked a mesh with nz of 1
            
            >>> nx = 1
            >>> ny = 2
            >>> nz = 3
            
            >>> mesh = UniformGrid3D(nx=nx, ny=ny, nz=nz, dx=dx, dy=dy, dz=dz)

            >>> cellVertexIDs = numerix.array((9, 8, 7, 6, 3, 2, 1, 0))
            >>> cellVertexIDs = numerix.array((cellVertexIDs, cellVertexIDs + 2, cellVertexIDs + 6,
            ...                                cellVertexIDs + 8, cellVertexIDs + 12, cellVertexIDs + 14))
            >>> cellVertexIDs = cellVertexIDs.swapaxes(0,1)
            >>> print parallel.procID > 0 or numerix.allclose(mesh._cellVertexIDs, cellVertexIDs)
            1
            
            >>> nx = 3
            >>> ny = 1
            >>> nz = 2
            
            >>> mesh = UniformGrid3D(nx=nx, ny=ny, nz=nz, dx=dx, dy=dy, dz=dz)

            >>> cellVertexIDs = numerix.array((13, 12, 9, 8, 5, 4, 1, 0))
            >>> cellVertexIDs = numerix.array((cellVertexIDs, cellVertexIDs + 1, cellVertexIDs + 2,
            ...                                cellVertexIDs + 8, cellVertexIDs + 9, cellVertexIDs + 10))
            >>> cellVertexIDs = cellVertexIDs.swapaxes(0,1)
            >>> print parallel.procID > 0 or numerix.allclose(mesh._cellVertexIDs, cellVertexIDs)
            1

        """

def _test():
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
