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

from fipy.tools import numerix
from fipy.tools.numerix import MA
from fipy.tools.decorators import getsetDeprecated
from fipy.tools import parallelComm

from fipy.meshes.uniformGrid import UniformGrid
from fipy.meshes.builders import _UniformGrid3DBuilder
from fipy.meshes.builders import _Grid3DBuilder
from fipy.meshes.representations.gridRepresentation import _Grid3DRepresentation
from fipy.meshes.topologies.gridTopology import _Grid3DTopology

__all__ = ["UniformGrid3D"]

class UniformGrid3D(UniformGrid):
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
                 origin = [[0], [0], [0]], overlap=2, communicator=parallelComm,
                 _RepresentationClass=_Grid3DRepresentation,
                 _TopologyClass=_Grid3DTopology):

        super(UniformGrid3D, self).__init__(communicator=communicator,
                                            _RepresentationClass=_RepresentationClass,
                                            _TopologyClass=_TopologyClass)

        builder = _UniformGrid3DBuilder()

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
        
        builder.buildGridData([dx, dy, dz], [nx, ny, nz], overlap, 
                              communicator, origin)
                                                        
        ([self.dx, self.dy, self.dz],
         [self.nx, self.ny, self.nz],
         self.dim,
         scale,
         self.globalNumberOfCells,
         self.globalNumberOfFaces,
         self.overlap,
         self.offset,
         self.numberOfVertices,
         self.numberOfFaces,
         self.numberOfCells,
         self.shape,
         self.physicalShape,
         self._meshSpacing,
         self.numberOfXYFaces,
         self.numberOfXZFaces,
         self.numberOfYZFaces,
         self.numberOfHorizontalRows,
         self.numberOfVerticalColumns,
         self.numberOfLayers,
         self.origin) = builder.gridData
        
    """
    Topology set and calc
    """

    @property
    def _exteriorFaces(self):
        """
        Return only the faces that have one neighboring cell.
        """
        XYids = self._XYFaceIDs
        XZids = self._XZFaceIDs
        YZids = self._YZFaceIDs
        
        exteriorIDs = numerix.concatenate((numerix.ravel(XYids[...,      0].swapaxes(0,1)), 
                                           numerix.ravel(XYids[...,     -1].swapaxes(0,1)),
                                           numerix.ravel(XZids[...,  0,...]), 
                                           numerix.ravel(XZids[..., -1,...]),
                                           numerix.ravel(YZids[ 0,     ...]), 
                                           numerix.ravel(YZids[-1,     ...])))
                                                     
        from fipy.variables.faceVariable import FaceVariable
        exteriorFaces = FaceVariable(mesh=self, value=False)
        exteriorFaces[exteriorIDs] = True
        return exteriorFaces

    @property
    def _interiorFaces(self):
        """
        Return only the faces that have two neighboring cells
        """
        XYids = self._XYFaceIDs
        XZids = self._XZFaceIDs
        YZids = self._YZFaceIDs
        
        interiorIDs = numerix.concatenate((numerix.ravel(XYids[ ...     ,1:-1]),
                                           numerix.ravel(XZids[ ...,1:-1, ...]),
                                           numerix.ravel(YZids[1:-1,      ...].swapaxes(0,1))))
                                                     
        from fipy.variables.faceVariable import FaceVariable
        interiorFaces = FaceVariable(mesh=self, value=False)
        interiorFaces[interiorIDs] = True
        return interiorFaces

    @property
    def _cellToFaceOrientations(self):
        tmp = numerix.take(self.faceCellIDs[0], self.cellFaceIDs)
        return (tmp == MA.indices(tmp.shape)[-1]) * 2 - 1

    @property
    def _adjacentCellIDs(self):
        faceCellIDs = self.faceCellIDs
        return (MA.where(MA.getmaskarray(faceCellIDs[0]), faceCellIDs[1], faceCellIDs[0]).filled(),
                MA.where(MA.getmaskarray(faceCellIDs[1]), faceCellIDs[0], faceCellIDs[1]).filled())

    @property
    def _cellToCellIDs(self):
        ids = MA.zeros((6, self.nx, self.ny, self.nz), 'l')
        indices = numerix.indices((self.nx, self.ny, self.nz))
        ids[0] = indices[0] + (indices[1] + indices[2] * self.ny) * self.nx - 1
        ids[1] = indices[0] + (indices[1] + indices[2] * self.ny) * self.nx + 1
        ids[2] = indices[0] + (indices[1] + indices[2] * self.ny - self.nz) * self.nx
        ids[3] = indices[0] + (indices[1] + indices[2] * self.ny + self.nz) * self.nx
        ids[4] = indices[0] + (indices[1] + (indices[2] - 1) * self.ny) * self.nx
        ids[5] = indices[0] + (indices[1] + (indices[2] + 1) * self.ny) * self.nx
        
        ids[0, 0,    ...] = MA.masked
        ids[1,-1,    ...] = MA.masked
        ids[2,..., 0,...] = MA.masked
        ids[3,...,-1,...] = MA.masked
        ids[4,...,     0] = MA.masked
        ids[5,...,    -1] = MA.masked

        return MA.reshape(ids.swapaxes(1,3), (6, self.numberOfCells))
        
    @property
    def _cellToCellIDsFilled(self):
        N = self.numberOfCells
        M = self._maxFacesPerCell
        cellIDs = numerix.repeat(numerix.arange(N)[numerix.newaxis, ...], M, axis=0)
        cellToCellIDs = self._cellToCellIDs
        return MA.where(MA.getmaskarray(cellToCellIDs), cellIDs, cellToCellIDs)     
                                                                                                
    """
    Geometry set and calc
    """

    @property
    def _faceAreas(self):
        return numerix.concatenate((numerix.repeat((self.dx * self.dy,), self.numberOfXYFaces),
                                    numerix.repeat((self.dx * self.dz,), self.numberOfXZFaces),
                                    numerix.repeat((self.dy * self.dz,), self.numberOfYZFaces)))

    @property
    def faceNormals(self):
        XYnor = numerix.zeros((3, self.nx, self.ny, self.nz + 1), 'l')
        XYnor[0,      ...] =  1
        XYnor[0,  ...,  0] = -1

        XZnor = numerix.zeros((3, self.nx, self.ny + 1, self.nz), 'l')
        XZnor[1,      ...] =  1
        XZnor[1,...,0,...] = -1

        YZnor = numerix.zeros((3, self.nx + 1, self.ny, self.nz), 'l')
        YZnor[2,      ...] =  1
        YZnor[2, 0,   ...] = -1
        
        return numerix.concatenate((numerix.reshape(XYnor[::-1].swapaxes(1,3), (3, self.numberOfXYFaces)), 
                                    numerix.reshape(XZnor[::-1].swapaxes(1,3), (3, self.numberOfXZFaces)), 
                                    numerix.reshape(YZnor[::-1].swapaxes(1,3), (3, self.numberOfYZFaces))), axis=1)

    @property
    def _cellVolumes(self):
        return numerix.ones(self.numberOfCells, 'd') * self.dx * self.dy * self.dz

    @property
    def _cellCenters(self):
        centers = numerix.zeros((3, self.nx, self.ny, self.nz), 'd')
        indices = numerix.indices((self.nx, self.ny, self.nz))
        centers[0] = (indices[0] + 0.5) * self.dx
        centers[1] = (indices[1] + 0.5) * self.dy
        centers[2] = (indices[2] + 0.5) * self.dz
        ccs = numerix.reshape(centers.swapaxes(1,3), (3, self.numberOfCells)) + self.origin
        return ccs

    @property
    def _cellDistances(self):
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

    @property
    def _faceToCellDistanceRatio(self):
        XYdis = numerix.zeros((self.nx, self.ny, self.nz + 1),'d')
        XYdis[:] = 0.5
        XYdis[..., 0] = 1
        XYdis[...,-1] = 1
        
        XZdis = numerix.zeros((self.nx, self.ny + 1, self.nz),'d')
        XZdis[:] = 0.5
        XZdis[..., 0,...] = 1
        XZdis[...,-1,...] = 1
        
        YZdis = numerix.zeros((self.nx + 1, self.ny, self.nz),'d')
        YZdis[:] = 0.5
        YZdis[ 0,...] = 1
        YZdis[-1,...] = 1
        
        return numerix.concatenate((numerix.ravel(XYdis.swapaxes(0,2)),
                                    numerix.ravel(XZdis.swapaxes(0,2)),
                                    numerix.ravel(YZdis.swapaxes(0,2))), axis=1)
    
    @property
    def _orientedFaceNormals(self):
        return self.faceNormals

    @property
    def _faceTangents1(self):
        XYtan = numerix.zeros((3, self.nx, self.ny, self.nz + 1), 'l')
        XYtan[2,      ...] =  1
        
        XZtan = numerix.zeros((3, self.nx, self.ny + 1, self.nz), 'l')
        XZtan[2,      ...] =  1
        
        YZtan = numerix.zeros((3, self.nx + 1, self.ny, self.nz), 'l')
        YZtan[1,      ...] =  1
        
        return numerix.concatenate((numerix.reshape(XYtan[::-1].swapaxes(1,3), (3, self.numberOfXYFaces)), 
                                    numerix.reshape(XZtan[::-1].swapaxes(1,3), (3, self.numberOfXZFaces)), 
                                    numerix.reshape(YZtan[::-1].swapaxes(1,3), (3, self.numberOfYZFaces))), axis=1)
        
    @property
    def _faceTangents2(self):
        XYtan = numerix.zeros((3, self.nx, self.ny, self.nz + 1), 'l')
        XYtan[1,      ...] =  1
        
        XZtan = numerix.zeros((3, self.nx, self.ny + 1, self.nz), 'l')
        XZtan[0,      ...] =  1
        
        YZtan = numerix.zeros((3, self.nx + 1, self.ny, self.nz), 'l')
        YZtan[0,      ...] =  1
        
        return numerix.concatenate((numerix.reshape(XYtan[::-1].swapaxes(1,3), (3, self.numberOfXYFaces)), 
                                    numerix.reshape(XZtan[::-1].swapaxes(1,3), (3, self.numberOfXZFaces)), 
                                    numerix.reshape(YZtan[::-1].swapaxes(1,3), (3, self.numberOfYZFaces))), axis=1)
    
    @property
    def _cellToCellDistances(self):
        distances = numerix.zeros((6, self.nx, self.ny, self.nz), 'd')
        distances[0] = self.dx
        distances[1] = self.dx
        distances[2] = self.dy
        distances[3] = self.dy
        distances[4] = self.dz
        distances[5] = self.dz
        
        distances[0,  0,...    ] = self.dx / 2.
        distances[1, -1,...    ] = self.dx / 2.
        distances[2,...,  0,...] = self.dy / 2.
        distances[3,..., -1,...] = self.dy / 2.
        distances[4,...,      0] = self.dz / 2.
        distances[5,...,     -1] = self.dz / 2.

        return numerix.reshape(distances.swapaxes(1,3), (self.numberOfCells, 6))
        
    @property
    def _cellNormals(self):
        normals = numerix.zeros((3, 6, self.numberOfCells), 'd')
        normals[...,0,...] = [[-1], [ 0], [ 0]]
        normals[...,1,...] = [[ 1], [ 0], [ 0]]
        normals[...,2,...] = [[ 0], [-1], [ 0]]
        normals[...,3,...] = [[ 0], [ 1], [ 0]]
        normals[...,4,...] = [[ 0], [ 0], [-1]]
        normals[...,5,...] = [[ 0], [ 0], [ 1]]

        return normals
        
    @property
    def _cellAreas(self):
        areas = numerix.ones((6, self.numberOfCells), 'd')
        areas[0] = self.dy * self.dz
        areas[1] = self.dy * self.dz
        areas[2] = self.dx * self.dz
        areas[3] = self.dx * self.dz
        areas[4] = self.dx * self.dy
        areas[5] = self.dx * self.dy
        return areas

    @property
    def _cellAreaProjections(self):
        return self._cellAreas * self._cellNormals

##         from numMesh/mesh

    @property
    def _faceCenters(self):
                                  
        XYcen = numerix.zeros((3, self.nx, self.ny, self.nz + 1), 'd')
        indices = numerix.indices((self.nx, self.ny, self.nz + 1))
        XYcen[0] = (indices[0] + 0.5) * self.dx
        XYcen[1] = (indices[1] + 0.5) * self.dy
        XYcen[2] = indices[2] * self.dz

        XZcen = numerix.zeros((3, self.nx, self.ny + 1, self.nz), 'd')
        indices = numerix.indices((self.nx, self.ny + 1, self.nz))
        XZcen[0] = (indices[0] + 0.5) * self.dx
        XZcen[1] = indices[1] * self.dy
        XZcen[2] = (indices[2] + 0.5) * self.dz
        
        YZcen = numerix.zeros((3, self.nx + 1, self.ny, self.nz), 'd')
        indices = numerix.indices((self.nx + 1, self.ny, self.nz))
        YZcen[0] = indices[0] * self.dx
        YZcen[1] = (indices[1] + 0.5) * self.dy
        YZcen[2] = (indices[2] + 0.5) * self.dz

        return numerix.concatenate((numerix.reshape(XYcen.swapaxes(1,3), (3, self.numberOfXYFaces)), 
                                    numerix.reshape(XZcen.swapaxes(1,3), (3, self.numberOfXZFaces)),
                                    numerix.reshape(YZcen.swapaxes(1,3), (3, self.numberOfYZFaces))), axis=1) + self.origin
                                                                 
    @property
    def _orientedAreaProjections(self):
        return self._areaProjections

    @property
    def _areaProjections(self):
        return self.faceNormals * self._faceAreas
     
    @property
    def _faceAspectRatios(self):
        return self._faceAreas / self._cellDistances  
         
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
        from fipy.meshes.nonUniformGrid3D import NonUniformGrid3D
        args = self.args.copy()
        origin = args['origin']
        from fipy.tools import serialComm
        args['communicator'] = serialComm
        del args['origin']
        return NonUniformGrid3D(**args) + origin

    @property
    def _cellFaceIDs(self):
        return MA.array(_Grid3DBuilder.createCells(self.nx,
                                                   self.ny,
                                                   self.nz,
                                                   self.numberOfXYFaces,
                                                   self.numberOfXZFaces,
                                                   self.numberOfYZFaces))

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

    @property
    def vertexCoords(self):
        return _Grid3DBuilder.createVertices(self.dx, self.dy, self.dz,
                                             self.nx, self.ny, self.nz,
                                             self.numberOfVertices,     
                                             self.numberOfHorizontalRows,
                                             self.numberOfVerticalColumns) \
                + self.origin

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
     
    @property
    def _cellVertexIDs(self):
        return self._orderedCellVertexIDs

    @property
    def faceVertexIDs(self):
       return _Grid3DBuilder.createFaces(self.nx, self.ny, self.nz)[1]

    def _calcOrderedCellVertexIDs(self):
        """Correct ordering for VTK_VOXEL"""
        ids = numerix.zeros((8, self.nx, self.ny, self.nz), 'l')
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
        
##     scaling
    
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
            >>> print numerix.allequal(XYFaceIDs, mesh._XYFaceIDs) # doctest: +PROCESSOR_0
            True
              
            >>> XZFaceIDs = numerix.array((((12,), (15,), (18,)),
            ...                            ((13,), (16,), (19,)),
            ...                            ((14,), (17,), (20,))))
            >>> print numerix.allequal(mesh._XZFaceIDs, XZFaceIDs) # doctest: +PROCESSOR_0
            True
            
            >>> YZFaceIDs = numerix.array((((21,), (25,)),
            ...                            ((22,), (26,)),
            ...                            ((23,), (27,)),
            ...                            ((24,), (28,))))
            >>> print numerix.allequal(mesh._YZFaceIDs, YZFaceIDs) # doctest: +PROCESSOR_0
            True

            >>> adjacentCellIDs = (numerix.array([0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 0, 1, 2, 3, 4, 5, 0, 0, 1, 2, 3, 3, 4, 5]), 
            ...                    numerix.array([0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 3, 4, 5, 0, 1, 2, 2, 3, 4, 5, 5]))
            >>> print numerix.allequal(mesh._adjacentCellIDs, adjacentCellIDs) # doctest: +PROCESSOR_0
            True

            >>> vertices = numerix.array(((0., 1., 2., 3., 0., 1., 2., 3., 0., 1., 2., 3., 0., 1., 2., 3., 0., 1., 2., 3., 0., 1., 2., 3.),
            ...                           (0., 0., 0., 0., 1., 1., 1., 1., 2., 2., 2., 2., 0., 0., 0., 0., 1., 1., 1., 1., 2., 2., 2., 2.),
            ...                           (0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.)))
            >>> vertices *= numerix.array([[dx], [dy], [dz]])
            
            >>> print numerix.allequal(vertices, mesh.vertexCoords) # doctest: +PROCESSOR_0
            True
        
            >>> faces = numerix.array((( 0,  1,  2,  4,  5,  6, 12, 13, 14, 16, 17, 18,  0,  1,  2,  4,  5,  6,  8,  9, 10,  0,  1,  2,  3,  4,  5,  6,  7),
            ...                        ( 1,  2,  3,  5,  6,  7, 13, 14, 15, 17, 18, 19,  1,  2,  3,  5,  6,  7,  9, 10, 11,  4,  5,  6,  7,  8,  9, 10, 11),
            ...                        ( 5,  6,  7,  9, 10, 11, 17, 18, 19, 21, 22, 23, 13, 14, 15, 17, 18, 19, 21, 22, 23, 16, 17, 18, 19, 20, 21, 22, 23),
            ...                        ( 4,  5,  6,  8,  9, 10, 16, 17, 18, 20, 21, 22, 12, 13, 14, 16, 17, 18, 20, 21, 22, 12, 13, 14, 15, 16, 17, 18, 19))) 
            >>> print numerix.allclose(faces, mesh.faceVertexIDs) # doctest: +PROCESSOR_0
            True

            >>> cells = numerix.array(((21, 22, 23, 25, 26, 27),
            ...                        (22, 23, 24, 26, 27, 28),
            ...                        (12, 13, 14, 15, 16, 17),
            ...                        (15, 16, 17, 18, 19, 20),
            ...                        ( 0,  1,  2,  3,  4,  5),
            ...                        ( 6,  7,  8,  9, 10, 11)))
            >>> print numerix.allequal(cells, mesh.cellFaceIDs) # doctest: +PROCESSOR_0
            True

            >>> externalFaces = numerix.array((0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 18, 19, 20, 21, 24, 25, 28))
            >>> print numerix.allequal(externalFaces, 
            ...                        numerix.nonzero(mesh.exteriorFaces)) # doctest: +PROCESSOR_0
            True

            >>> internalFaces = numerix.array((15, 16, 17, 22, 23, 26, 27))
            >>> print numerix.allequal(internalFaces, 
            ...                        numerix.nonzero(mesh.interiorFaces)) # doctest: +PROCESSOR_0
            True

            >>> from fipy.tools.numerix import MA
            >>> faceCellIds = MA.masked_values((( 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 0, 1, 2, 3, 4, 5, 0, 0, 1, 2, 3, 3, 4, 5),
            ...                                 (-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 3, 4, 5,-1,-1,-1,-1, 1, 2,-1,-1, 4, 5,-1)), -1) 
            >>> print numerix.allequal(faceCellIds, mesh.faceCellIDs) # doctest: +PROCESSOR_0
            True
            
            >>> xy = dx * dy
            >>> xz = dx * dz
            >>> yz = dy * dz
            >>> faceAreas = numerix.array((xy, xy, xy, xy, xy, xy, xy, xy, xy, xy, xy, xy,
            ...                            xz, xz, xz, xz, xz, xz, xz, xz, xz,
            ...                            yz, yz, yz, yz, yz, yz, yz, yz))
            >>> print numerix.allclose(faceAreas, mesh._faceAreas, atol = 1e-10, rtol = 1e-10) # doctest: +PROCESSOR_0
            1
            
            >>> faceCoords = numerix.take(vertices, faces, axis=1)
            >>> faceCenters = (faceCoords[...,0,:] + faceCoords[...,1,:] + faceCoords[...,2,:] + faceCoords[...,3,:]) / 4.
            >>> print numerix.allclose(faceCenters, mesh.faceCenters, atol = 1e-10, rtol = 1e-10)
            True

            >>> faceNormals = numerix.array((( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 1, 1,-1, 1, 1, 1),
            ...                              ( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0),
            ...                              (-1,-1,-1,-1,-1,-1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)))
            >>> print numerix.allclose(faceNormals, mesh.faceNormals, atol = 1e-10, rtol = 1e-10) # doctest: +PROCESSOR_0
            True

            >>> cellToFaceOrientations = numerix.array((( 1,-1,-1, 1,-1,-1),
            ...                                         ( 1, 1, 1, 1, 1, 1),
            ...                                         ( 1, 1, 1,-1,-1,-1),
            ...                                         ( 1, 1, 1, 1, 1, 1),
            ...                                         ( 1, 1, 1, 1, 1, 1),
            ...                                         ( 1, 1, 1, 1, 1, 1)))
            >>> print numerix.allequal(cellToFaceOrientations, mesh._cellToFaceOrientations) # doctest: +PROCESSOR_0
            True
                                             
            >>> cellVolumes = numerix.array((dx*dy*dz, dx*dy*dz, dx*dy*dz, dx*dy*dz, dx*dy*dz, dx*dy*dz))
            >>> print numerix.allclose(cellVolumes, mesh.cellVolumes, atol = 1e-10, rtol = 1e-10) # doctest: +PROCESSOR_0
            True

            >>> cellCenters = numerix.array(((   dx/2., 3.*dx/2., 5.*dx/2.,   dx/2., 3.*dx/2., 5.*dx/2.),
            ...                              (   dy/2.,    dy/2.,    dy/2.,3.*dy/2., 3.*dy/2., 3.*dy/2.),
            ...                              (   dz/2.,    dz/2.,    dz/2.,   dz/2.,    dz/2.,    dz/2.)))
            >>> print numerix.allclose(cellCenters, mesh.cellCenters, atol = 1e-10, rtol = 1e-10)
            True
                                              
            >>> cellDistances = numerix.array((dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2,
            ...                                dy/2, dy/2, dy/2, dy, dy, dy, dy/2, dy/2, dy/2,
            ...                                dx/2, dx, dx, dx/2, dx/2, dx, dx, dx/2))
            >>> print numerix.allclose(cellDistances, mesh._cellDistances, atol = 1e-10, rtol = 1e-10) # doctest: +PROCESSOR_0
            True
            
            >>> faceToCellDistances = MA.masked_values(((dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dy/2, dy/2, dy/2, dy/2, dy/2, dy/2, dy/2, dy/2, dy/2, dx/2, dx/2, dx/2, dx/2, dx/2, dx/2, dx/2, dx/2),
            ...                                         (  -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1, dy/2, dy/2, dy/2,   -1,   -1,   -1,   -1, dx/2, dx/2,   -1,   -1, dx/2, dx/2,   -1)), -1)
            >>> faceToCellDistanceRatios = faceToCellDistances[0] / cellDistances
            >>> print numerix.allclose(faceToCellDistanceRatios, mesh._faceToCellDistanceRatio, atol = 1e-10, rtol = 1e-10) # doctest: +PROCESSOR_0
            True

            >>> areaProjections = faceNormals * faceAreas
            >>> print numerix.allclose(areaProjections, mesh._areaProjections, atol = 1e-10, rtol = 1e-10) # doctest: +PROCESSOR_0
            True

            >>> tangents1 = numerix.array(((1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0),
            ...                            (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1),
            ...                            (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)))
            >>> print numerix.allclose(tangents1, mesh._faceTangents1, atol = 1e-10, rtol = 1e-10) # doctest: +PROCESSOR_0
            True

            >>> tangents2 = numerix.array(((0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
            ...                            (1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
            ...                            (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)))
            >>> print numerix.allclose(tangents2, mesh._faceTangents2, atol = 1e-10, rtol = 1e-10) # doctest: +PROCESSOR_0
            True

            >>> cellToCellIDs = MA.masked_values(((-1, 0, 1, -1, 3, 4),
            ...                                   (1, 2, -1, 4, 5, -1),
            ...                                   (-1, -1, -1, 0, 1, 2),
            ...                                   (3, 4, 5, -1, -1, -1),
            ...                                   (-1, -1, -1, -1, -1, -1),
            ...                                   (-1, -1, -1, -1, -1, -1)), -1)
            >>> print numerix.allequal(mesh._cellToCellIDs, cellToCellIDs) # doctest: +PROCESSOR_0
            True

            >>> cellToCellIDsFilled = numerix.array([[0, 0, 1, 3, 3, 4],
            ...                                      [1, 2, 2, 4, 5, 5],
            ...                                      [0, 1, 2, 0, 1, 2],
            ...                                      [3, 4, 5, 3, 4, 5],
            ...                                      [0, 1, 2, 3, 4, 5],
            ...                                      [0, 1, 2, 3, 4, 5]])
            >>> print numerix.allequal(mesh._cellToCellIDsFilled, cellToCellIDsFilled) # doctest: +PROCESSOR_0
            True
              
            >>> cellToCellDistances = numerix.take(cellDistances, cells)
            >>> print numerix.allclose(cellToCellDistances, mesh._cellToCellDistances, atol = 1e-10, rtol = 1e-10) # doctest: +PROCESSOR_0
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
            >>> print numerix.allclose(cellNormals, mesh._cellNormals, atol = 1e-10, rtol = 1e-10) # doctest: +PROCESSOR_0
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
            >>> print numerix.allclose(cellAreaProjections, mesh._cellAreaProjections, atol = 1e-10, rtol = 1e-10) # doctest: +PROCESSOR_0
            True

            >>> cellVertexIDs = numerix.array((17, 16, 13, 12, 5, 4, 1, 0))
            >>> cellVertexIDs = numerix.array((cellVertexIDs, cellVertexIDs + 1, cellVertexIDs + 2,
            ...                                cellVertexIDs + 4, cellVertexIDs + 5, cellVertexIDs + 6))
            >>> cellVertexIDs = cellVertexIDs.swapaxes(0,1)
            >>> print numerix.allclose(mesh._cellVertexIDs, cellVertexIDs) # doctest: +PROCESSOR_0
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
            >>> print numerix.allclose(mesh._cellVertexIDs, cellVertexIDs) # doctest: +PROCESSOR_0
            1
            
            >>> nx = 3
            >>> ny = 1
            >>> nz = 2
            
            >>> mesh = UniformGrid3D(nx=nx, ny=ny, nz=nz, dx=dx, dy=dy, dz=dz)

            >>> cellVertexIDs = numerix.array((13, 12, 9, 8, 5, 4, 1, 0))
            >>> cellVertexIDs = numerix.array((cellVertexIDs, cellVertexIDs + 1, cellVertexIDs + 2,
            ...                                cellVertexIDs + 8, cellVertexIDs + 9, cellVertexIDs + 10))
            >>> cellVertexIDs = cellVertexIDs.swapaxes(0,1)
            >>> print numerix.allclose(mesh._cellVertexIDs, cellVertexIDs) # doctest: +PROCESSOR_0
            1

        """

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
