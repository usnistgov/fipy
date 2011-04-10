#!/usr/bin/env python

## 
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "mesh.py"
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
 #  See the file "license.terms" for information on usage and  redistribution
 #  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #  
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

from fipy.tools import numerix
from fipy.tools.numerix import MA

from fipy.variables.faceVariable import FaceVariable
from fipy.variables.cellVariable import CellVariable

from abstractGeometries import AbstractScaledMeshGeometry
from abstractUniformGeometries import AbstractUniformGridGeometry

class UniformGridScaledGeometry3D(AbstractScaledMeshGeometry):

    def __init__(self, geom):
        self._geom = geom
                                   
    @property
    def orientedAreaProjections(self):
        return self.areaProjections

    @property
    def areaProjections(self):
        return self._geom.faceNormals * self._geom.faceAreas
     
    @property
    def faceAspectRatios(self):
        return self._geom.faceAreas / self._geom.cellDistances
    
class UniformGridGeometry3D(AbstractUniformGridGeometry):
    
    def __init__(self, mesh,
                       dx, dy, dz,
                       nx, ny, nz,
                       numberOfCells, 
                       numberOfXYFaces, numberOfXZFaces, numberOfYZFaces,
                       origin):

        self.mesh = mesh
        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.numberOfCells = numberOfCells
        self.numberOfXYFaces = numberOfXYFaces
        self.numberOfXZFaces = numberOfXZFaces
        self.numberOfYZFaces = numberOfYZFaces
        self.origin = origin

        self._scaledGeometry = UniformGridScaledGeometry3D(self)
                                 
    @property
    def faceAreas(self):
        return FaceVariable(mesh=self.mesh, 
                            value=numerix.concatenate((numerix.repeat((self.dx * self.dy,), self.numberOfXYFaces),
                                                       numerix.repeat((self.dx * self.dz,), self.numberOfXZFaces),
                                                       numerix.repeat((self.dy * self.dz,), self.numberOfYZFaces))))

    @property
    def faceNormals(self):
        XYnor = numerix.zeros((3, self.nx, self.ny, self.nz + 1))
        XYnor[0,      ...] =  1
        XYnor[0,  ...,  0] = -1

        XZnor = numerix.zeros((3, self.nx, self.ny + 1, self.nz))
        XZnor[1,      ...] =  1
        XZnor[1,...,0,...] = -1

        YZnor = numerix.zeros((3, self.nx + 1, self.ny, self.nz))
        YZnor[2,      ...] =  1
        YZnor[2, 0,   ...] = -1
        
        return FaceVariable(mesh=self.mesh,
                            value=numerix.concatenate((numerix.reshape(XYnor[::-1].swapaxes(1,3), (3, self.numberOfXYFaces)), 
                                                       numerix.reshape(XZnor[::-1].swapaxes(1,3), (3, self.numberOfXZFaces)), 
                                                       numerix.reshape(YZnor[::-1].swapaxes(1,3), (3, self.numberOfYZFaces))), axis=1))

    @property
    def faceCellToCellNormals(self):
        return self.faceNormals
        
    @property
    def cellVolumes(self):
        return CellVariable(mesh=self.mesh, value=self.dx * self.dy * self.dz)

    @property
    def cellCenters(self):
        centers = numerix.zeros((3, self.nx, self.ny, self.nz), 'd')
        indices = numerix.indices((self.nx, self.ny, self.nz))
        centers[0] = (indices[0] + 0.5) * self.dx
        centers[1] = (indices[1] + 0.5) * self.dy
        centers[2] = (indices[2] + 0.5) * self.dz
        return CellVariable(mesh=self.mesh,
                            value=numerix.reshape(centers.swapaxes(1,3), (3, self.numberOfCells)) + self.origin)

    @property
    def cellDistances(self):
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

        return FaceVariable(mesh=self.mesh, 
                            value=numerix.concatenate((numerix.ravel(XYdis),
                                                       numerix.ravel(XZdis),
                                                       numerix.ravel(YZdis))))

    @property
    def faceToCellDistanceRatio(self):
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
        
        return FaceVariable(mesh=self.mesh, 
                            value=numerix.concatenate((numerix.ravel(XYdis.swapaxes(0,2)),
                                                       numerix.ravel(XZdis.swapaxes(0,2)),
                                                       numerix.ravel(YZdis.swapaxes(0,2))), axis=1))
    
    @property
    def orientedFaceNormals(self):
        return self.faceNormals

    @property
    def faceTangents1(self):
        XYtan = numerix.zeros((3, self.nx, self.ny, self.nz + 1))
        XYtan[2,      ...] =  1
        
        XZtan = numerix.zeros((3, self.nx, self.ny + 1, self.nz))
        XZtan[2,      ...] =  1
        
        YZtan = numerix.zeros((3, self.nx + 1, self.ny, self.nz))
        YZtan[1,      ...] =  1
        
        return FaceVariable(mesh=self.mesh, 
                            value=numerix.concatenate((numerix.reshape(XYtan[::-1].swapaxes(1,3), (3, self.numberOfXYFaces)), 
                                                       numerix.reshape(XZtan[::-1].swapaxes(1,3), (3, self.numberOfXZFaces)), 
                                                       numerix.reshape(YZtan[::-1].swapaxes(1,3), (3, self.numberOfYZFaces))), axis=1))
        
    @property
    def faceTangents2(self):
        XYtan = numerix.zeros((3, self.nx, self.ny, self.nz + 1))
        XYtan[1,      ...] =  1
        
        XZtan = numerix.zeros((3, self.nx, self.ny + 1, self.nz))
        XZtan[0,      ...] =  1
        
        YZtan = numerix.zeros((3, self.nx + 1, self.ny, self.nz))
        YZtan[0,      ...] =  1
        
        return FaceVariable(mesh=self.mesh, 
                            value=numerix.concatenate((numerix.reshape(XYtan[::-1].swapaxes(1,3), (3, self.numberOfXYFaces)), 
                                                       numerix.reshape(XZtan[::-1].swapaxes(1,3), (3, self.numberOfXZFaces)), 
                                                       numerix.reshape(YZtan[::-1].swapaxes(1,3), (3, self.numberOfYZFaces))), axis=1))
    
    @property
    def cellToCellDistances(self):
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

        return CellVariable(mesh=self.mesh, 
                            value=numerix.reshape(distances.swapaxes(1,3), (self.numberOfCells, 6)),
                            elementshape=(6,))
        
    @property
    def cellNormals(self):
        normals = CellVariable(mesh=self.mesh, value=0., elementshape=(3,6))
        normals[...,0,...] = [[-1], [ 0], [ 0]]
        normals[...,1,...] = [[ 1], [ 0], [ 0]]
        normals[...,2,...] = [[ 0], [-1], [ 0]]
        normals[...,3,...] = [[ 0], [ 1], [ 0]]
        normals[...,4,...] = [[ 0], [ 0], [-1]]
        normals[...,5,...] = [[ 0], [ 0], [ 1]]

        return normals
        
    @property
    def cellAreas(self):
        areas = CellVariable(mesh=self.mesh, value=0., elementshape=(6,))
        areas[0] = self.dy * self.dz
        areas[1] = self.dy * self.dz
        areas[2] = self.dx * self.dz
        areas[3] = self.dx * self.dz
        areas[4] = self.dx * self.dy
        areas[5] = self.dx * self.dy
        return areas

    @property
    def cellAreaProjections(self):
        return self.cellAreas * self.cellNormals

##         from numMesh/mesh

    @property
    def faceCenters(self):
                                  
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

        return FaceVariable(mesh=self.mesh,
                            value=(numerix.concatenate((numerix.reshape(XYcen.swapaxes(1,3), (3, self.numberOfXYFaces)), 
                                                       numerix.reshape(XZcen.swapaxes(1,3), (3, self.numberOfXZFaces)),
                                                       numerix.reshape(YZcen.swapaxes(1,3), (3, self.numberOfYZFaces))), axis=1)
                                   + self.origin),
                            rank=1)
     
