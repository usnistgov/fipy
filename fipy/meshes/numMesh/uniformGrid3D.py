#!/usr/bin/env python

## -*-Pyth-*-
 # ########################################################################
 # FiPy - a finite volume PDE solver in Python
 # 
 # FILE: "uniformGrid3D.py"
 #                                     created: 3/2/06 {3:57:15 PM}
 #                                 last update: 3/7/06 {4:59:52 PM}
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
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
 # History
 # 
 # modified   by  rev reason
 # ---------- --- --- -----------
 # 2006-03-02 JEG 1.0 original
 # 
 # ########################################################################
 ##

import MA

from fipy.meshes.numMesh.grid3D import Grid3D
from fipy.meshes.meshIterator import FaceIterator
from fipy.tools import numerix
from fipy.tools.dimensions.physicalField import PhysicalField

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
    def __init__(self, dx = 1., dy = 1., dz = 1., nx = 1, ny = 1, nz = 1, origin = (0,0,0)):
        self.dim = 3
        
        self.dx = PhysicalField(value = dx)
        scale = PhysicalField(value = 1, unit = self.dx.getUnit())
        self.dx /= scale
        
        self.nx = nx
        
        self.dy = PhysicalField(value = dy)
        if self.dy.getUnit().isDimensionless():
            self.dy = dy
        else:
            self.dy /= scale
            
        self.ny = ny
        
        self.dz = PhysicalField(value = dy)
        if self.dz.getUnit().isDimensionless():
            self.dz = dz
        else:
            self.dz /= scale
            
        self.nz = nz
        
        self.origin = PhysicalField(value = origin)
        self.origin /= scale

        self.numberOfVertices = (self.nx + 1) * (self.ny + 1) * (self.nz + 1)
        self.numberOfXYFaces = self.nx * self.ny * (self.nz + 1)
        self.numberOfXZFaces = self.nx * (self.ny + 1) * self.nz
        self.numberOfYZFaces = (self.nx + 1) * self.ny * self.nz
        self.numberOfFaces = self.numberOfXYFaces + self.numberOfXZFaces + self.numberOfYZFaces
        self.numberOfCells = self.nx * self.ny * self.nz
        
        
        self.scale = {
            'length': 1.,
            'area': 1.,
            'volume': 1.
        }

        self.setScale(value = scale)
        
    def _translate(self, vector):
        return UniformGrid3D(dx = self.dx, nx = self.nx, 
                             dy = self.dy, ny = self.ny, 
                             dz = self.dz, nz = self.nz, 
                             origin = self.origin + vector)

    def __mul__(self, factor):
        return UniformGrid3D(dx = self.dx * factor, nx = self.nx, 
                             dy = self.dy * factor, ny = self.ny, 
                             dz = self.dz * factor, nz = self.nz, 
                             origin = self.origin * factor)

    def _getConcatenableMesh(self):
        from fipy.meshes.numMesh.mesh3D import Mesh3D
        return Mesh3D(vertexCoords = self.getVertexCoords(), 
                      faceVertexIDs = self._createFaces(), 
                      cellFaceIDs = self._createCells())
                      
    def _concatenate(self, other, smallNumber):
        return self._getConcatenableMesh()._concatenate(other = other, smallNumber = smallNumber)
        
##     get topology methods

##         from common/mesh
        
    def _getCellFaceIDs(self):
        return MA.array(self._createCells())
        
    def _getXYFaceIDs(self):
        ids = numerix.arange(0, self.numberOfXYFaces)
        return numerix.reshape(ids, (self.nz + 1, self.ny, self.nx))

    def _getXZFaceIDs(self):
        ids = numerix.arange(self.numberOfXYFaces, self.numberOfXYFaces + self.numberOfXZFaces)
        return numerix.reshape(ids, (self.nz, self.ny + 1, self.nx))

    def _getYZFaceIDs(self):
        ids = numerix.arange(self.numberOfXYFaces + self.numberOfXZFaces, self.numberOfFaces)
        return numerix.reshape(ids, (self.nz, self.ny, self.nx + 1))

    def getExteriorFaces(self):
        XYids = self._getXYFaceIDs()
        XZids = self._getXZFaceIDs()
        YZids = self._getYZFaceIDs()
        return FaceIterator(mesh=self,
                            ids=numerix.concatenate((numerix.ravel(XYids[  0,    ...]), 
                                                     numerix.ravel(XYids[ -1,    ...]),
                                                     numerix.ravel(XZids[...,  0,...]), 
                                                     numerix.ravel(XZids[..., -1,...]),
                                                     numerix.ravel(YZids[...    ,  0]), 
                                                     numerix.ravel(YZids[...    , -1]))))
        
    def getInteriorFaces(self):
        XYids = self._getXYFaceIDs()
        XZids = self._getXZFaceIDs()
        YZids = self._getYZFaceIDs()
        return FaceIterator(mesh=self,
                            ids=numerix.concatenate((numerix.ravel(XYids[1:-1,      ...]),
                                                     numerix.ravel(XZids[ ...,1:-1, ...]),
                                                     numerix.ravel(YZids[ ...     ,1:-1]))))

    def _getCellFaceOrientations(self):
        tmp = numerix.MAtake(self.getFaceCellIDs()[...,0], self._getCellFaceIDs())
        return (tmp == MA.indices(tmp.shape)[0]) * 2 - 1

    def _getAdjacentCellIDs(self):
        faceCellIDs = self.getFaceCellIDs()
        return (MA.where(MA.getmaskarray(faceCellIDs[...,0]), faceCellIDs[...,1], faceCellIDs[...,0]).filled(),
                MA.where(MA.getmaskarray(faceCellIDs[...,1]), faceCellIDs[...,0], faceCellIDs[...,1]).filled())

    def _getCellToCellIDs(self):
        ids = MA.zeros((self.nz, self.ny, self.nx, 6))
        indices = numerix.indices((self.nz, self.ny, self.nx))
        ids[...,0] = indices[2] + (indices[1] + indices[0] * self.ny) * self.nx - 1
        ids[...,1] = indices[2] + (indices[1] + indices[0] * self.ny) * self.nx + 1
        ids[...,2] = indices[2] + (indices[1] + indices[0] * self.ny - self.nz) * self.nx
        ids[...,3] = indices[2] + (indices[1] + indices[0] * self.ny + self.nz) * self.nx
        ids[...,4] = indices[2] + (indices[1] + (indices[0] - 1) * self.ny) * self.nx
        ids[...,5] = indices[2] + (indices[1] + (indices[0] + 1) * self.ny) * self.nx
        
        ids[...,     0,0] = MA.masked
        ids[...,    -1,1] = MA.masked
        ids[..., 0,...,2] = MA.masked
        ids[...,-1,...,3] = MA.masked
        ids[ 0,    ...,4] = MA.masked
        ids[-1,    ...,5] = MA.masked

        return MA.reshape(ids, (self.numberOfCells, 6))
        
    def _getCellToCellIDsFilled(self):
        N = self.getNumberOfCells()
        M = self._getMaxFacesPerCell()
        cellIDs = numerix.reshape(numerix.repeat(numerix.arange(N), M), (N, M))
        cellToCellIDs = self._getCellToCellIDs()
        return MA.where(MA.getmaskarray(cellToCellIDs), cellIDs, cellToCellIDs)
        
    def _getMaxFacesPerCell(self):
        return 6
        
##         from numMesh/mesh

    def getVertexCoords(self):
        return self._createVertices() + self.origin

    def getFaceCellIDs(self):
        XYids = MA.zeros((self.nz + 1, self.ny, self.nx, 2))
        indices = numerix.indices((self.nz + 1, self.ny, self.nx))
        XYids[...,1] = indices[2] + (indices[1] + indices[0] * self.ny) * self.nx
        XYids[...,0] = XYids[...,1] - self.nx * self.ny
        XYids[ 0,...,0] = XYids[ 0,...,1]
        XYids[ 0,...,1] = MA.masked
        XYids[-1,...,1] = MA.masked
        
        XZids = MA.zeros((self.nz, self.ny + 1, self.nx, 2))
        indices = numerix.indices((self.nz, self.ny + 1, self.nx))
        XZids[...,1] = indices[2] + (indices[1] + indices[0] * self.ny) * self.nx
        XZids[...,0] = XZids[...,1] - self.nx
        XZids[..., 0,...,0] = XZids[..., 0,...,1]
        XZids[..., 0,...,1] = MA.masked
        XZids[...,-1,...,1] = MA.masked

        YZids = MA.zeros((self.nz, self.ny, self.nx + 1, 2))
        indices = numerix.indices((self.nz, self.ny, self.nx + 1))
        YZids[...,1] = indices[2] + (indices[1] + indices[0] * self.ny) * self.nx
        YZids[...,0] = YZids[...,1] - 1
        YZids[..., 0,0] = YZids[..., 0,1]
        YZids[..., 0,1] = MA.masked
        YZids[...,-1,1] = MA.masked

        return MA.concatenate((MA.reshape(XYids, (self.numberOfXYFaces, 2)), 
                               MA.reshape(XZids, (self.numberOfXZFaces, 2)), 
                               MA.reshape(YZids, (self.numberOfYZFaces, 2))))

##     get geometry methods
        
##         from common/mesh
        
    def _getFaceAreas(self):
        return numerix.concatenate((numerix.repeat((self.dx * self.dy,), self.numberOfXYFaces),
                                    numerix.repeat((self.dx * self.dz,), self.numberOfXZFaces),
                                    numerix.repeat((self.dy * self.dz,), self.numberOfYZFaces)))

    def _getFaceNormals(self):
        XYnor = numerix.zeros((self.nz + 1, self.ny, self.nx, 3))
        XYnor[      ...,0] =  1
        XYnor[0,    ...,0] = -1

        XZnor = numerix.zeros((self.nz, self.ny + 1, self.nx, 3))
        XZnor[      ...,1] =  1
        XZnor[...,0,...,1] = -1

        YZnor = numerix.zeros((self.nz, self.ny, self.nx + 1, 3))
        YZnor[      ...,2] =  1
        YZnor[    ...,0,2] = -1
        
        return numerix.concatenate((numerix.reshape(XYnor[...,::-1], (self.numberOfXYFaces, 3)), 
                                    numerix.reshape(XZnor[...,::-1], (self.numberOfXZFaces, 3)), 
                                    numerix.reshape(YZnor[...,::-1], (self.numberOfYZFaces, 3))))
        
    def getCellVolumes(self):
        return numerix.ones(self.numberOfCells, 'd') * self.dx * self.dy * self.dz

    def getCellCenters(self):
        centers = numerix.zeros((self.nz, self.ny, self.nx, 3), 'd')
        indices = numerix.indices((self.nz, self.ny, self.nx))
        centers[...,2] = (indices[2] + 0.5) * self.dx
        centers[...,1] = (indices[1] + 0.5) * self.dy
        centers[...,0] = (indices[0] + 0.5) * self.dz
        return numerix.reshape(centers[...,::-1], (self.numberOfCells, 3)) + self.origin

    def _getCellDistances(self):
        XYdis = numerix.zeros((self.nz + 1, self.ny, self.nx),'d')
        XYdis[:] = self.dz
        XYdis[ 0,...] = self.dz / 2.
        XYdis[-1,...] = self.dz / 2.
        
        XZdis = numerix.zeros((self.nz, self.ny + 1, self.nx),'d')
        XZdis[:] = self.dy
        XZdis[..., 0,...] = self.dy / 2.
        XZdis[...,-1,...] = self.dy / 2.

        YZdis = numerix.zeros((self.nz, self.ny, self.nx + 1),'d')
        YZdis[:] = self.dx
        YZdis[..., 0] = self.dx / 2.
        YZdis[...,-1] = self.dx / 2.

        return numerix.concatenate((numerix.ravel(XYdis),
                                    numerix.ravel(XZdis),
                                    numerix.ravel(YZdis)))

    def _getFaceToCellDistanceRatio(self):
        XYdis = numerix.zeros((self.nz + 1, self.ny, self.nx),'d')
        XYdis[:] = 0.5
        XYdis[ 0,...] = 1
        XYdis[-1,...] = 1
        
        XZdis = numerix.zeros((self.nz, self.ny + 1, self.nx),'d')
        XZdis[:] = 0.5
        XZdis[..., 0,...] = 1
        XZdis[...,-1,...] = 1
        
        YZdis = numerix.zeros((self.nz, self.ny, self.nx + 1),'d')
        YZdis[:] = 0.5
        YZdis[..., 0] = 1
        YZdis[...,-1] = 1
        
        return numerix.concatenate((numerix.ravel(XYdis),
                                    numerix.ravel(XZdis),
                                    numerix.ravel(YZdis)))
                                    
    def _getOrientedAreaProjections(self):
        return self._getAreaProjections()

    def _getAreaProjections(self):
        return self._getFaceNormals() * self._getFaceAreas()[..., numerix.NewAxis]

    def _getOrientedFaceNormals(self):
        return self._getFaceNormals()

    def _getFaceTangents1(self):
        XYtan = numerix.zeros((self.nz + 1, self.ny, self.nx, 3))
        XYtan[      ...,2] =  1
        
        XZtan = numerix.zeros((self.nz, self.ny + 1, self.nx, 3))
        XZtan[      ...,2] =  1
        
        YZtan = numerix.zeros((self.nz, self.ny, self.nx + 1, 3))
        YZtan[      ...,1] =  1
        
        return numerix.concatenate((numerix.reshape(XYtan[...,::-1], (self.numberOfXYFaces, 3)), 
                                    numerix.reshape(XZtan[...,::-1], (self.numberOfXZFaces, 3)), 
                                    numerix.reshape(YZtan[...,::-1], (self.numberOfYZFaces, 3))))
        
    def _getFaceTangents2(self):
        XYtan = numerix.zeros((self.nz + 1, self.ny, self.nx, 3))
        XYtan[      ...,1] =  1
        
        XZtan = numerix.zeros((self.nz, self.ny + 1, self.nx, 3))
        XZtan[      ...,0] =  1
        
        YZtan = numerix.zeros((self.nz, self.ny, self.nx + 1, 3))
        YZtan[      ...,0] =  1
        
        return numerix.concatenate((numerix.reshape(XYtan[...,::-1], (self.numberOfXYFaces, 3)), 
                                    numerix.reshape(XZtan[...,::-1], (self.numberOfXZFaces, 3)), 
                                    numerix.reshape(YZtan[...,::-1], (self.numberOfYZFaces, 3))))
        
    def _getFaceAspectRatios(self):
        return self._getFaceAreas() / self._getCellDistances()
    
    def _getCellToCellDistances(self):
        distances = numerix.zeros((self.nz, self.ny, self.nx, 6), 'd')
        distances[...,0] = self.dx
        distances[...,1] = self.dx
        distances[...,2] = self.dy
        distances[...,3] = self.dy
        distances[...,4] = self.dz
        distances[...,5] = self.dz
        
        distances[...,      0,0] = self.dx / 2.
        distances[...,     -1,1] = self.dx / 2.
        distances[...,  0,...,2] = self.dy / 2.
        distances[..., -1,...,3] = self.dy / 2.
        distances[  0,...    ,4] = self.dz / 2.
        distances[ -1,...    ,5] = self.dz / 2.

        return numerix.reshape(distances, (self.numberOfCells, 6))
        
    def _getCellNormals(self):
        normals = numerix.zeros((self.numberOfCells, 6, 3), 'd')
        normals[...,0,...] = (-1, 0, 0)
        normals[...,1,...] = ( 1, 0, 0)
        normals[...,2,...] = ( 0,-1, 0)
        normals[...,3,...] = ( 0, 1, 0)
        normals[...,4,...] = ( 0, 0,-1)
        normals[...,5,...] = ( 0, 0, 1)

        return normals
        
    def _getCellAreas(self):
        areas = numerix.ones((self.numberOfCells,6), 'd')
        areas[...,0] = self.dy * self.dz
        areas[...,1] = self.dy * self.dz
        areas[...,2] = self.dx * self.dz
        areas[...,3] = self.dx * self.dz
        areas[...,4] = self.dx * self.dy
        areas[...,5] = self.dx * self.dy
        return areas

    def _getCellAreaProjections(self):
        return self._getCellAreas()[...,numerix.NewAxis] * self._getCellNormals()

##         from numMesh/mesh

    def getFaceCenters(self):
        XYcen = numerix.zeros((self.nz + 1, self.ny, self.nx, 3), 'd')
        indices = numerix.indices((self.nz + 1, self.ny, self.nx))
        XYcen[...,0] = indices[0] * self.dz
        XYcen[...,1] = (indices[1] + 0.5) * self.dy
        XYcen[...,2] = (indices[2] + 0.5) * self.dx

        XZcen = numerix.zeros((self.nz, self.ny + 1, self.nx, 3), 'd')
        indices = numerix.indices((self.nz, self.ny + 1, self.nx))
        XZcen[...,0] = (indices[0] + 0.5) * self.dz
        XZcen[...,1] = indices[1] * self.dy
        XZcen[...,2] = (indices[2] + 0.5) * self.dx
        
        YZcen = numerix.zeros((self.nz, self.ny, self.nx + 1, 3), 'd')
        indices = numerix.indices((self.nz, self.ny, self.nx + 1))
        YZcen[...,0] = (indices[0] + 0.5) * self.dz
        YZcen[...,1] = (indices[1] + 0.5) * self.dy
        YZcen[...,2] = indices[2] * self.dx

        return numerix.concatenate((numerix.reshape(XYcen[...,::-1], (self.numberOfXYFaces, 3)), 
                                    numerix.reshape(XZcen[...,::-1], (self.numberOfXZFaces, 3)),
                                    numerix.reshape(YZcen[...,::-1], (self.numberOfYZFaces, 3))))
                                    
    def _getCellVertexIDs(self):
        ids = numerix.zeros((self.nz, self.ny, self.nx, 8))
        indices = numerix.indices((self.nz, self.ny, self.nx))
        ids[...,1] = indices[2] + (indices[1] + (indices[0] + 1) * (self.ny + 1) + 1) * (self.nx + 1)
        ids[...,0] = ids[...,1] + 1
        ids[...,3] = indices[2] + (indices[1] + (indices[0] + 1) * (self.ny + 1)) * (self.nx + 1)
        ids[...,2] = ids[...,3] + 1
        ids[...,5] = indices[2] + (indices[1] + indices[0] * self.ny + 1) * (self.nx + 1)
        ids[...,4] = ids[...,5] + 1
        ids[...,7] = indices[2] + (indices[1] + indices[0] * self.ny) * (self.nx + 1)
        ids[...,6] = ids[...,7] + 1
        
        return numerix.reshape(ids, (self.numberOfCells, 8))
        
    def _getFaceVertexIDs(self):
        return self._createFaces()
                                    
    def _getOrderedCellVertexIDs(self):
        """Correct ordering for VTK_VOXEL"""
        return self._getCellVertexIDs()     
        
##     scaling
    
    def _calcScaledGeometry(self):
        pass
        
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
            
            >>> mesh = UniformGrid3D(nx = nx, ny = ny, nz = nz, dx = dx, dy = dy, dz = dz)
            
            >>> print mesh._getXYFaceIDs()
            [[[ 0, 1, 2,]
              [ 3, 4, 5,]]
             [[ 6, 7, 8,]
              [ 9,10,11,]]]
              
            >>> print mesh._getXZFaceIDs()
            [ [[12,13,14,]
              [15,16,17,]
              [18,19,20,]]]
              
            >>> print mesh._getYZFaceIDs()
            [ [[21,22,23,24,]
              [25,26,27,28,]]]

            >>> print mesh._getAdjacentCellIDs()
            ([0,1,2,3,4,5,0,1,2,3,4,5,0,1,2,0,1,2,3,4,5,0,0,1,2,3,3,4,5,], [0,1,2,3,4,5,0,1,2,3,4,5,0,1,2,3,4,5,3,4,5,0,1,2,2,3,4,5,5,])

            >>> vertices = numerix.array(((0., 0., 0.), (1., 0., 0.), (2., 0., 0.), (3., 0., 0.),
            ...                           (0., 1., 0.), (1., 1., 0.), (2., 1., 0.), (3., 1., 0.),
            ...                           (0., 2., 0.), (1., 2., 0.), (2., 2., 0.), (3., 2., 0.),
            ...                           (0., 0., 1.), (1., 0., 1.), (2., 0., 1.), (3., 0., 1.),
            ...                           (0., 1., 1.), (1., 1., 1.), (2., 1., 1.), (3., 1., 1.),
            ...                           (0., 2., 1.), (1., 2., 1.), (2., 2., 1.), (3., 2., 1.)))
            >>> vertices *= numerix.array((dx, dy, dz))
            
            >>> numerix.allequal(vertices, mesh._createVertices())
            1
        
            >>> faces = numerix.array(((0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (4, 5, 9, 8), (5, 6, 10, 9), (6, 7, 11, 10),
            ...                        (12, 13, 17, 16), (13, 14, 18, 17), (14, 15, 19, 18), (16, 17, 21, 20), (17, 18, 22, 21), (18, 19, 23, 22),
            ...                        (0, 1, 13, 12), (1, 2, 14, 13), (2, 3, 15, 14), (4, 5, 17, 16), (5, 6, 18, 17), (6, 7, 19, 18), (8, 9, 21, 20), (9, 10, 22, 21), (10, 11, 23, 22),
            ...                        (0, 4, 16, 12), (1, 5, 17, 13), (2, 6, 18, 14), (3, 7, 19, 15), (4, 8, 20, 16), (5, 9, 21, 17), (6, 10, 22, 18), (7, 11, 23, 19)))
            >>> numerix.allequal(faces, mesh._createFaces())
            1

            >>> cells = numerix.array(((21, 22, 12, 15, 0, 6),
            ...                       (22, 23, 13, 16, 1, 7),
            ...                       (23, 24, 14, 17, 2, 8),
            ...                       (25, 26, 15, 18, 3, 9),
            ...                       (26, 27, 16, 19, 4, 10),
            ...                       (27, 28, 17, 20, 5, 11)))
            >>> numerix.allequal(cells, mesh._createCells())
            1

            >>> externalFaces = numerix.array((0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 18, 19, 20, 21, 25, 24, 28))
            >>> numerix.allequal(externalFaces, mesh.getExteriorFaces())
            1

            >>> internalFaces = numerix.array((15, 16, 17, 22, 23, 26, 27))
            >>> numerix.allequal(internalFaces, mesh.getInteriorFaces())
            1

            >>> import MA
            >>> faceCellIds = MA.masked_values(((0, -1), (1, -1), (2, -1), (3, -1), (4, -1), (5, -1),
            ...                                 (0, -1), (1, -1), (2, -1), (3, -1), (4, -1), (5, -1),
            ...                                 (0, -1), (1, -1), (2, -1), (0, 3), (1, 4), (2, 5), (3, -1), (4, -1), (5, -1),
            ...                                 (0, -1), (0, 1), (1, 2), (2, -1), (3, -1), (3, 4), (4, 5), (5, -1)), -1)
            >>> numerix.allequal(faceCellIds, mesh.getFaceCellIDs())
            1
            
            >>> xy = dx * dy
            >>> xz = dx * dz
            >>> yz = dy * dz
            >>> faceAreas = numerix.array((xy, xy, xy, xy, xy, xy, xy, xy, xy, xy, xy, xy,
            ...                            xz, xz, xz, xz, xz, xz, xz, xz, xz,
            ...                            yz, yz, yz, yz, yz, yz, yz, yz))
            >>> numerix.allclose(faceAreas, mesh._getFaceAreas(), atol = 1e-10, rtol = 1e-10)
            1
            
            >>> faceCoords = numerix.take(vertices, faces)
            >>> faceCenters = (faceCoords[:,0] + faceCoords[:,1] + faceCoords[:,2] + faceCoords[:, 3]) / 4.
            >>> numerix.allclose(faceCenters, mesh.getFaceCenters(), atol = 1e-10, rtol = 1e-10)
            1

            >>> faceNormals = numerix.array(((0, 0, -1), (0, 0, -1), (0, 0, -1), (0, 0, -1), (0, 0, -1), (0, 0, -1),
            ...                              (0, 0, 1), (0, 0, 1), (0, 0, 1), (0, 0, 1), (0, 0, 1), (0, 0, 1),
            ...                              (0, -1, 0), (0, -1, 0), (0, -1, 0), (0, 1, 0), (0, 1, 0), (0, 1, 0), (0, 1, 0), (0, 1, 0), (0, 1, 0),
            ...                              (-1, 0, 0), (1, 0, 0), (1, 0, 0), (1, 0, 0), (-1, 0, 0), (1, 0, 0), (1, 0, 0), (1, 0, 0)))
            >>> numerix.allclose(faceNormals, mesh._getFaceNormals(), atol = 1e-10, rtol = 1e-10)
            1

            >>> cellToFaceOrientations = numerix.array(((1, 1, 1, 1, 1, 1),
            ...                                         (-1, 1, 1, 1, 1, 1),
            ...                                         (-1, 1, 1, 1, 1, 1),
            ...                                         (1, 1, -1, 1, 1, 1),
            ...                                         (-1, 1, -1, 1, 1, 1),
            ...                                         (-1, 1, -1, 1, 1, 1)))
            >>> numerix.allequal(cellToFaceOrientations, mesh._getCellFaceOrientations())
            1
                                             
            >>> cellVolumes = numerix.array((dx*dy*dz, dx*dy*dz, dx*dy*dz, dx*dy*dz, dx*dy*dz, dx*dy*dz))
            >>> numerix.allclose(cellVolumes, mesh.getCellVolumes(), atol = 1e-10, rtol = 1e-10)
            1

            >>> cellCenters = numerix.array(((dx/2.,dy/2.,dz/2.), (3.*dx/2.,dy/2.,dz/2.), (5.*dx/2.,dy/2.,dz/2.),
            ...                              (dx/2.,3.*dy/2.,dz/2.), (3.*dx/2.,3.*dy/2.,dz/2.), (5.*dx/2.,3.*dy/2.,dz/2.)))
            >>> numerix.allclose(cellCenters, mesh.getCellCenters(), atol = 1e-10, rtol = 1e-10)
            1
                                              
            >>> cellDistances = numerix.array((dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2,
            ...                                dy/2, dy/2, dy/2, dy, dy, dy, dy/2, dy/2, dy/2,
            ...                                dx/2, dx, dx, dx/2, dx/2, dx, dx, dx/2))
            >>> numerix.allclose(cellDistances, mesh._getCellDistances(), atol = 1e-10, rtol = 1e-10)
            1
            
            >>> faceToCellDistances = MA.masked_values(((dz/2, -1), (dz/2, -1), (dz/2, -1), (dz/2, -1), (dz/2, -1), (dz/2, -1),
            ...                                         (dz/2, -1), (dz/2, -1), (dz/2, -1), (dz/2, -1), (dz/2, -1), (dz/2, -1),
            ...                                         (dy/2, -1), (dy/2, -1), (dy/2, -1), (dy/2, dy/2), (dy/2, dy/2), (dy/2, dy/2), (dy/2, -1), (dy/2, -1), (dy/2, -1),
            ...                                         (dx/2, -1), (dx/2, dx/2), (dx/2, dx/2), (dx/2, -1), (dx/2, -1), (dx/2, dx/2), (dx/2, dx/2), (dx/2, -1)), -1)
            >>> faceToCellDistanceRatios = faceToCellDistances[...,0] / cellDistances
            >>> numerix.allclose(faceToCellDistanceRatios, mesh._getFaceToCellDistanceRatio(), atol = 1e-10, rtol = 1e-10)
            1

            >>> areaProjections = faceNormals * faceAreas[...,numerix.NewAxis]
            >>> numerix.allclose(areaProjections, mesh._getAreaProjections(), atol = 1e-10, rtol = 1e-10)
            1

            >>> tangents1 = numerix.array(((1, 0, 0), (1, 0, 0), (1, 0, 0), (1, 0, 0), (1, 0, 0), (1, 0, 0), (1, 0, 0), (1, 0, 0), (1, 0, 0), (1, 0, 0), (1, 0, 0), (1, 0, 0),
            ...                            (1, 0, 0), (1, 0, 0), (1, 0, 0), (1, 0, 0), (1, 0, 0), (1, 0, 0), (1, 0, 0), (1, 0, 0), (1, 0, 0),
            ...                            (0, 1, 0), (0, 1, 0), (0, 1, 0), (0, 1, 0), (0, 1, 0), (0, 1, 0), (0, 1, 0), (0, 1, 0)))
            >>> numerix.allclose(tangents1, mesh._getFaceTangents1(), atol = 1e-10, rtol = 1e-10)
            1

            >>> tangents2 = numerix.array(((0, 1, 0), (0, 1, 0), (0, 1, 0), (0, 1, 0), (0, 1, 0), (0, 1, 0), (0, 1, 0), (0, 1, 0), (0, 1, 0), (0, 1, 0), (0, 1, 0), (0, 1, 0),
            ...                            (0, 0, 1), (0, 0, 1), (0, 0, 1), (0, 0, 1), (0, 0, 1), (0, 0, 1), (0, 0, 1), (0, 0, 1), (0, 0, 1),
            ...                            (0, 0, 1), (0, 0, 1), (0, 0, 1), (0, 0, 1), (0, 0, 1), (0, 0, 1), (0, 0, 1), (0, 0, 1)))
            >>> numerix.allclose(tangents2, mesh._getFaceTangents2(), atol = 1e-10, rtol = 1e-10)
            1

            >>> print mesh._getCellToCellIDs()
            [[-- ,1 ,-- ,3 ,-- ,-- ,]
             [0 ,2 ,-- ,4 ,-- ,-- ,]
             [1 ,-- ,-- ,5 ,-- ,-- ,]
             [-- ,4 ,0 ,-- ,-- ,-- ,]
             [3 ,5 ,1 ,-- ,-- ,-- ,]
             [4 ,-- ,2 ,-- ,-- ,-- ,]]

            >>> print mesh._getCellToCellIDsFilled()
            [[0,1,0,3,0,0,]
             [0,2,1,4,1,1,]
             [1,2,2,5,2,2,]
             [3,4,0,3,3,3,]
             [3,5,1,4,4,4,]
             [4,5,2,5,5,5,]]
              
            >>> cellToCellDistances = numerix.take(cellDistances, cells)
            >>> numerix.allclose(cellToCellDistances, mesh._getCellToCellDistances(), atol = 1e-10, rtol = 1e-10)
            1

            >>> cellNormals = numerix.array((((-1, 0, 0), (1, 0, 0), (0, -1, 0), (0, 1, 0), (0, 0, -1), (0, 0, 1)),
            ...                              ((-1, 0, 0), (1, 0, 0), (0, -1, 0), (0, 1, 0), (0, 0, -1), (0, 0, 1)),
            ...                              ((-1, 0, 0), (1, 0, 0), (0, -1, 0), (0, 1, 0), (0, 0, -1), (0, 0, 1)),
            ...                              ((-1, 0, 0), (1, 0, 0), (0, -1, 0), (0, 1, 0), (0, 0, -1), (0, 0, 1)),
            ...                              ((-1, 0, 0), (1, 0, 0), (0, -1, 0), (0, 1, 0), (0, 0, -1), (0, 0, 1)),
            ...                              ((-1, 0, 0), (1, 0, 0), (0, -1, 0), (0, 1, 0), (0, 0, -1), (0, 0, 1)) ))
            >>> numerix.allclose(cellNormals, mesh._getCellNormals(), atol = 1e-10, rtol = 1e-10)
            1

            >>> vv = numerix.array(((-yz, 0, 0), (yz, 0, 0), (0, -xz, 0), (0, xz, 0), (0, 0, -xy), (0, 0, xy)))
            >>> cellAreaProjections = numerix.array(((vv,vv,vv,vv,vv,vv)))
            >>> numerix.allclose(cellAreaProjections, mesh._getCellAreaProjections(), atol = 1e-10, rtol = 1e-10)
            1

            >>> cellVertexIDs = numerix.array((17, 16, 13, 12, 5, 4, 1, 0))
            >>> cellVertexIDs = numerix.array((cellVertexIDs, cellVertexIDs + 1, cellVertexIDs + 2,
            ...                                cellVertexIDs + 4, cellVertexIDs + 5, cellVertexIDs + 6))
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
