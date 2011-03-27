#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "grid3D.py"
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

__docformat__ = 'restructuredtext'

from fipy.tools import numerix

from fipy.meshes.mesh import Mesh
from fipy.tools import vector
from fipy.tools.dimensions.physicalField import PhysicalField
from fipy.tools.decorators import getsetDeprecated

from fipy.tools import parallel

from fipy.meshes.builders import NonuniformGrid3DBuilder
from fipy.meshes.gridlike import Gridlike3D

class Grid3D(Mesh):
    """
    3D rectangular-prism Mesh

    X axis runs from left to right.
    Y axis runs from bottom to top.
    Z axis runs from front to back.

    Numbering System:

    Vertices: Numbered in the usual way. X coordinate changes most quickly, then Y, then Z.

    Cells: Same numbering system as vertices.

    Faces: XY faces numbered first, then XZ faces, then YZ faces. Within each subcategory, it is numbered in the usual way.
    """
    def __init__(self, dx = 1., dy = 1., dz = 1., nx = None, ny = None, nz = None, overlap=2, communicator=parallel):

        builder = NonuniformGrid3DBuilder()
        
        self.args = {
            'dx': dx, 
            'dy': dy,
            'dz' :dz,
            'nx': nx, 
            'ny': ny,
            'nz': nz,
            'overlap': overlap,
            'communicator': communicator
        }
        
        builder.buildGridData([dx, dy, dz], [nx, ny, nz], overlap, 
                              communicator)
                                                                      
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
         self.numberOfLayersDeep,
         vertices,
         faces,
         cells,
         self.Xoffset, self.Yoffset, self.Zoffset) = builder.gridData
        
        Mesh.__init__(self, vertices, faces, cells)
        
        self._setScale(scaleLength = scale)
         
    def __getstate__(self):
        return Gridlike3D.__getstate__(self)

    def __setstate__(self, dict):
        return Gridlike3D.__setstate__(self, dict)

    def __repr__(self):
        return Gridlike3D.__repr__(self)

    def _isOrthogonal(self):
        return Gridlike3D._isOrthogonal(self)

    @property
    def _concatenatedClass(self):
        return Gridlike3D._concatenatedClass
                                                                
    @property
    def _globalNonOverlappingCellIDs(self):
        """
        Return the IDs of the local mesh in the context of the
        global parallel mesh. Does not include the IDs of boundary cells.

        E.g., would return [0, 1, 4, 5] for mesh A

            A        B
        ------------------
        | 4 | 5 || 6 | 7 |
        ------------------
        | 0 | 1 || 2 | 3 |
        ------------------
        
        .. note:: Trivial except for parallel meshes
        """
        return Gridlike3D._globalNonOverlappingCellIDs(self)

    @property
    def _globalOverlappingCellIDs(self):
        """
        Return the IDs of the local mesh in the context of the
        global parallel mesh. Includes the IDs of boundary cells.
        
        E.g., would return [0, 1, 2, 4, 5, 6] for mesh A

            A        B
        ------------------
        | 4 | 5 || 6 | 7 |
        ------------------
        | 0 | 1 || 2 | 3 |
        ------------------
        
        .. note:: Trivial except for parallel meshes
        """
        return Gridlike3D._globalOverlappingCellIDs(self)

    @property
    def _localNonOverlappingCellIDs(self):
        """
        Return the IDs of the local mesh in isolation. 
        Does not include the IDs of boundary cells.
        
        E.g., would return [0, 1, 2, 3] for mesh A

            A        B
        ------------------
        | 3 | 4 || 4 | 5 |
        ------------------
        | 0 | 1 || 1 | 2 |
        ------------------
        
        .. note:: Trivial except for parallel meshes
        """
        return Gridlike3D._localNonOverlappingCellIDs(self)

    @property
    def _localOverlappingCellIDs(self):
        """
        Return the IDs of the local mesh in isolation. 
        Includes the IDs of boundary cells.
        
        E.g., would return [0, 1, 2, 3, 4, 5] for mesh A

            A        B
        ------------------
        | 3 | 4 || 5 |   |
        ------------------
        | 0 | 1 || 2 |   |
        ------------------
        
        .. note:: Trivial except for parallel meshes
        """
        return Gridlike3D._localOverlappingCellIDs(self)
 
    def _calcScaleArea(self):
        return self.scale['length']**2

    def _calcScaleVolume(self):
        return self.scale['length']**3  

    def _calcFaceNormals(self):
        XYFaceNormals = numerix.zeros((3, self.numberOfXYFaces))
        XYFaceNormals[2, (self.nx * self.ny):] = 1
        XYFaceNormals[2, :(self.nx * self.ny)] = -1
        XZFaceNormals = numerix.zeros((3, self.numberOfXZFaces))
        xzd = numerix.arange(self.numberOfXZFaces)
        xzd = xzd % (self.nx * (self.ny + 1))
        xzd = (xzd < self.nx)
        xzd = 1 - (2 * xzd)
        XZFaceNormals[1, :] = xzd
        YZFaceNormals = numerix.zeros((3, self.numberOfYZFaces))
        YZFaceNormals[0, :] = 1
        YZFaceNormals[0, ::self.nx + 1] = -1
        return numerix.concatenate((XYFaceNormals, 
                                    XZFaceNormals, 
                                    YZFaceNormals), 
                                   axis=-1)
        
    def _calcFaceTangents(self):
        ## need to see whether order matters.
        faceTangents1 = numerix.zeros((3, self.numberOfFaces), 'd')
        faceTangents2 = numerix.zeros((3, self.numberOfFaces), 'd')
        ## XY faces
        faceTangents1[0, :self.numberOfXYFaces] = 1.
        faceTangents2[1, :self.numberOfXYFaces] = 1.
        ## XZ faces
        faceTangents1[0, self.numberOfXYFaces:self.numberOfXYFaces + self.numberOfXZFaces] = 1.
        faceTangents2[2, self.numberOfXYFaces:self.numberOfXYFaces + self.numberOfXZFaces] = 1.
        ## YZ faces
        faceTangents1[1, self.numberOfXYFaces + self.numberOfXZFaces:] = 1.
        faceTangents2[2, self.numberOfXYFaces + self.numberOfXZFaces:] = 1.
        return faceTangents1, faceTangents2                                     

## The following method is broken when dx, dy or dz are not scalar. Simpler to use the generic
## _calcFaceAreas rather than do the required type checking, resizing and outer product.
##
##     def _calcFaceAreas(self):
##         XYFaceAreas = numerix.ones(self.numberOfXYFaces)
##         XYFaceAreas = XYFaceAreas * self.dx * self.dy
##         XZFaceAreas = numerix.ones(self.numberOfXZFaces)
##         XZFaceAreas = XZFaceAreas * self.dx * self.dz        
##         YZFaceAreas = numerix.ones(self.numberOfYZFaces)
##         YZFaceAreas = YZFaceAreas * self.dy * self.dz
##         self.faceAreas =  numerix.concatenate((XYFaceAreas, XZFaceAreas, YZFaceAreas))

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
            
            >>> mesh = Grid3D(nx = nx, ny = ny, nz = nz, dx = dx, dy = dy, dz = dz)

            >>> adjacentCellIDs = (numerix.array([0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 0,
            ...                               1, 2, 3, 4, 5, 0, 0, 1, 2, 3, 3, 4, 5]),
            ...                numerix.array([0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3,
            ...                               4, 5, 3, 4, 5, 0, 1, 2, 2, 3, 4, 5, 5]))
            >>> print parallel.procID > 0 or numerix.allequal(mesh._adjacentCellIDs, adjacentCellIDs)
            True

            >>> vertices = numerix.array(((0., 1., 2., 3., 0., 1., 2., 3., 0., 1., 2., 3., 
            ...                            0., 1., 2., 3., 0., 1., 2., 3., 0., 1., 2., 3.),
            ...                           (0., 0., 0., 0., 1., 1., 1., 1., 2., 2., 2., 2.,
            ...                            0., 0., 0., 0., 1., 1., 1., 1., 2., 2., 2., 2.),
            ...                           (0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
            ...                            1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.)))
            >>> vertices *= numerix.array([[dx], [dy], [dz]])
            
            >>> print parallel.procID > 0 or numerix.allequal(vertices,
            ...                                               mesh.vertexCoords)
            True
        
            >>> faces = numerix.array(((0, 1, 2, 4,  5,  6, 12, 13, 14, 16, 17, 18,  0,  1,  2,  4,  5,  6,  8,  9, 10,  0,  1,  2,  3,  4,  5,  6,  7),
            ...                        (1, 2, 3, 5,  6,  7, 13, 14, 15, 17, 18, 19,  1,  2,  3,  5,  6,  7,  9, 10, 11,  4,  5,  6,  7,  8,  9, 10, 11),
            ...                        (5, 6, 7, 9, 10, 11, 17, 18, 19, 21, 22, 23, 13, 14, 15, 17, 18, 19, 21, 22, 23, 16, 17, 18, 19, 20, 21, 22, 23),
            ...                        (4, 5, 6, 8,  9, 10, 16, 17, 18, 20, 21, 22, 12, 13, 14, 16, 17, 18, 20, 21, 22, 12, 13, 14, 15, 16, 17, 18, 19)))
            >>> print parallel.procID > 0 or numerix.allequal(faces,
            ...                                               mesh.faceVertexIDs)
            True

            >>> cells = numerix.array(((21, 22, 23, 25, 26, 27),
            ...                        (22, 23, 24, 26, 27, 28),
            ...                        (12, 13, 14, 15, 16, 17),
            ...                        (15, 16, 17, 18, 19, 20),
            ...                        ( 0,  1,  2,  3,  4,  5),
            ...                        ( 6,  7,  8,  9, 10, 11)))
            >>> print parallel.procID > 0 or numerix.allequal(cells,
            ...                                               mesh.cellFaceIDs)
            True

            >>> externalFaces = numerix.array((0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 18, 19, 20, 21, 24, 25, 28))
            >>> print parallel.procID > 0 or numerix.allequal(externalFaces, 
            ...                              numerix.nonzero(mesh.exteriorFaces))
            True

            >>> internalFaces = numerix.array((15, 16, 17, 22, 23, 26, 27))
            >>> print parallel.procID > 0 or  numerix.allequal(internalFaces, 
            ...                               numerix.nonzero(mesh.interiorFaces))
            True

            >>> from fipy.tools.numerix import MA
            >>> faceCellIds = MA.masked_values((( 0,  1,  2,  3,  4,  5,  0,  1,  2,  3,  4,  5,  0,  1,  2,  0, 1, 2,  3,  4,  5,  0, 0, 1,  2,  3, 3, 4,  5),
            ...                                 (-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  3, 4, 5, -1, -1, -1, -1, 1, 2, -1, -1, 4, 5, -1)), -1)
            >>> print parallel.procID > 0 or numerix.allequal(faceCellIds, mesh.faceCellIDs)
            True
            
            >>> xy = dx * dy
            >>> xz = dx * dz
            >>> yz = dy * dz
            >>> faceAreas = numerix.array((xy, xy, xy, xy, xy, xy, xy, xy, xy, xy, xy, xy,
            ...                            xz, xz, xz, xz, xz, xz, xz, xz, xz,
            ...                            yz, yz, yz, yz, yz, yz, yz, yz))
            >>> print parallel.procID > 0 or numerix.allclose(faceAreas, mesh._faceAreas, atol = 1e-10, rtol = 1e-10)
            True
            
            >>> faceCoords = numerix.take(vertices, faces, axis=1)
            >>> faceCenters = (faceCoords[...,0,:] + faceCoords[...,1,:] + faceCoords[...,2,:] + faceCoords[...,3,:]) / 4.
            >>> print parallel.procID > 0 or numerix.allclose(faceCenters, mesh.faceCenters, atol = 1e-10, rtol = 1e-10)
            True

            >>> faceNormals = numerix.array((( 0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0,  0,  0,  0, 0, 0, 0, 0, 0, 0, -1, 1, 1, 1, -1, 1, 1, 1),
            ...                              ( 0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, -1, -1, -1, 1, 1, 1, 1, 1, 1,  0, 0, 0, 0,  0, 0, 0, 0),
            ...                              (-1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1,  0,  0,  0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0)))
            >>> print parallel.procID > 0 or numerix.allclose(faceNormals, mesh._faceNormals, atol = 1e-10, rtol = 1e-10)
            True

            >>> cellToFaceOrientations = numerix.array(((1, -1, -1, 1, -1, -1),
            ...                                         (1, 1, 1, 1, 1, 1),
            ...                                         (1, 1, 1, -1, -1, -1),
            ...                                         (1, 1, 1, 1, 1, 1),
            ...                                         (1, 1, 1, 1, 1, 1),
            ...                                         (1, 1, 1, 1, 1, 1)))
            >>> print parallel.procID > 0 or numerix.allequal(cellToFaceOrientations, mesh._cellToFaceOrientations)
            True
                                             
            >>> cellVolumes = numerix.array((dx*dy*dz, dx*dy*dz, dx*dy*dz, dx*dy*dz, dx*dy*dz, dx*dy*dz))
            >>> print parallel.procID > 0 or numerix.allclose(cellVolumes, mesh.cellVolumes, atol = 1e-10, rtol = 1e-10)
            True

            >>> cellCenters = numerix.array(((dx/2., 3.*dx/2., 5.*dx/2.,    dx/2., 3.*dx/2., 5.*dx/2.),
            ...                              (dy/2.,    dy/2.,    dy/2., 3.*dy/2., 3.*dy/2., 3.*dy/2.),
            ...                              (dz/2.,    dz/2.,    dz/2.,    dz/2.,    dz/2.,    dz/2.)))
            >>> print numerix.allclose(cellCenters, mesh.cellCenters, atol = 1e-10, rtol = 1e-10)
            True
                                              
            >>> faceToCellDistances = MA.masked_values(((dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dy/2, dy/2, dy/2, dy/2, dy/2, dy/2, dy/2, dy/2, dy/2, dx/2, dx/2, dx/2, dx/2, dx/2, dx/2, dx/2, dx/2),
            ...                                         (  -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1, dy/2, dy/2, dy/2,   -1,   -1,   -1,   -1, dx/2, dx/2,   -1,   -1, dx/2,   -1)), -1) 
            >>> print parallel.procID > 0 or numerix.allclose(faceToCellDistances, mesh._faceToCellDistances, atol = 1e-10, rtol = 1e-10)
            1
                                              
            >>> cellDistances = numerix.array((dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2,
            ...                                dy/2, dy/2, dy/2, dy, dy, dy, dy/2, dy/2, dy/2,
            ...                                dx/2, dx, dx, dx/2, dx/2, dx, dx, dx/2))
            >>> print parallel.procID > 0 or numerix.allclose(cellDistances, mesh._cellDistances, atol = 1e-10, rtol = 1e-10)
            True
            
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

            >>> cellToCellIDsFilled = numerix.array(((0, 0, 1, 3, 3, 4),
            ...                                      (1, 2, 2, 4, 5, 5),
            ...                                      (0, 1, 2, 0, 1, 2),
            ...                                      (3, 4, 5, 3, 4, 5),
            ...                                      (0, 1, 2, 3, 4, 5),
            ...                                      (0, 1, 2, 3, 4, 5)))            
            >>> print parallel.procID > 0 or numerix.allequal(mesh._cellToCellIDsFilled, cellToCellIDsFilled)
            True

            >>> cellToCellDistances = numerix.take(cellDistances, cells)
            >>> print parallel.procID > 0 or numerix.allclose(cellToCellDistances, mesh._cellToCellDistances, atol = 1e-10, rtol = 1e-10)
            True

            >>> interiorCellIDs = numerix.array(())
            >>> numerix.allequal(interiorCellIDs, mesh._interiorCellIDs)
            True

            >>> exteriorCellIDs = numerix.array((0, 1, 2, 3, 4, 5))
            >>> print parallel.procID > 0 or numerix.allequal(exteriorCellIDs, mesh._exteriorCellIDs)
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

            >>> cellAreaProjections = numerix.array((((-yz, -yz, -yz, -yz, -yz, -yz),
            ...                                       ( yz,  yz,  yz,  yz,  yz,  yz),
            ...                                       (  0,   0,   0,   0,   0,   0),
            ...                                       (  0,   0,   0,   0,   0,   0),
            ...                                       (  0,   0,   0,   0,   0,   0),
            ...                                       (  0,   0,   0,   0,   0,   0)),
            ...                                      ((  0,   0,   0,   0,   0,   0),
            ...                                       (  0,   0,   0,   0,   0,   0),
            ...                                       (-xz, -xz, -xz, -xz, -xz, -xz),
            ...                                       ( xz,  xz,  xz,  xz,  xz,  xz),
            ...                                       (  0,   0,   0,   0,   0,   0),
            ...                                       (  0,   0,   0,   0,   0,   0)),
            ...                                      ((  0,   0,   0,   0,   0,   0),
            ...                                       (  0,   0,   0,   0,   0,   0),
            ...                                       (  0,   0,   0,   0,   0,   0),
            ...                                       (  0,   0,   0,   0,   0,   0),
            ...                                       (-xy, -xy, -xy, -xy, -xy, -xy),
            ...                                       ( xy,  xy,  xy,  xy,  xy,  xy))))
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

            >>> print numerix.allclose(mesh.cellCenters, unpickledMesh.cellCenters)
            True

            The following test was for a bug when dx, dy or dz are arrays.
            The _calcFaceAreas() method was commented out to fix this.

            >>> Grid3D(nx=2., ny=2., nz=2., dx=(1., 2.), dy=(1., 2.), dz=(1., 2.))
            Grid3D(dx=(1.0, 2.0), dy=(1.0, 2.0), dz=(1.0, 2.0), nx=2, ny=2, nz=2)
        """

def _test():
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
