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
 #  Author: Alexander Mont <alexander.mont@nist.gov>
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
from fipy.tools import inline

from fipy.meshes.geometries.abstractGeometries import AbstractMeshGeometry
from fipy.meshes.geometries.abstractGeometries import AbstractScaledMeshGeometry

class UniformMeshScaledGeometry1D(AbstractScaledMeshGeometry):

    def __init__(self, geom):
        self.geom = geom
    
    @property
    def faceToCellDistanceRatio(self):
        distances = numerix.ones(self.geom.numberOfFaces, 'd')
        distances *= 0.5
        if len(distances) > 0:
            distances[0] = 1
            distances[-1] = 1
        return distances

    @property
    def areaProjections(self):
        return self.geom.faceNormals
        
    @property
    def orientedAreaProjections(self):
        return self._getAreaProjections()
        
    @property
    def _getFaceAspectRatios(self):
        return 1. / self._getCellDistances()

class UniformMeshScaledGeometry2D(AbstractScaledMeshGeometry):

    def __init__(self, geom, numFaces, dx, dy, nx, numHorizFaces):
        self.geom = geom
        self.numberOfFaces = numFaces
        self.dx = dx
        self.dy = dy
        self.nx = nx
        self.numberOfHorizontalFaces = numHorizFaces
      
    @property
    def orientedAreaProjections(self):
        return self._getAreaProjections()

    @property
    def areaProjections(self):
        return inline._optionalInline(self._getAreaProjectionsIn, self._getAreaProjectionsPy)

    def _getAreaProjectionsPy(self):
        return self.geom.faceNormals * self.geom.faceAreas

    def _getAreaProjectionsIn(self):
        areaProjections = numerix.zeros((2, self.numberOfFaces), 'd')

        inline._runInline("""
            if (i < nx) {
                areaProjections[i + 1 * ni] = -dx;
            } else if (i < Nhor) {
                areaProjections[i + 1 * ni] = dx;
            } else if ( (i - Nhor) % (nx + 1) == 0 ) {
                areaProjections[i + 0 * ni] = -dy;
            } else {
                areaProjections[i + 0 * ni] = dy;
           }
        """,
        dx = float(self.dx), # horrible hack to get around
        dy = float(self.dy), # http://www.scipy.org/scipy/scipy/ticket/496
        nx = self.nx,
        Nhor = self.numberOfHorizontalFaces,
        areaProjections = areaProjections,
        ni = self.numberOfFaces)

        return areaProjections
        
    @property
    def faceAspectRatios(self):
        return self.geom.faceAreas / self.geom.cellDistances
      
class UniformMeshGeometry(AbstractMeshGeometry):

    """Wrapped scaled geometry properties"""
    @property
    def scaledFaceAreas(self):
        return self.faceAreas

    @property
    def scaledCellVolumes(self):
        return self.cellVolumes

    @property
    def scaledCellCenters(self):
        return self.cellCenters

    @property
    def scaledCellDistances(self):
        return self.cellDistances

    @property
    def scaledCellToCellDistances(self):
        return self.cellToCellDistances

    @property
    def scaledFaceToCellDistances(self):
        return self.faceToCellDistances

    @property
    def faceToCellDistanceRatio(self):
        return self._scaledGeometry.faceToCellDistanceRatio

    @property
    def areaProjections(self):
        return self._scaledGeometry.areaProjections

    @property
    def orientedAreaProjections(self):
        return self._scaledGeometry.orientedAreaProjections

    """Geometry properties common to 1D, 2D, 3D"""
    
    @property
    def orientedFaceNormals(self):
        return self.faceNormals
      
class UniformMeshGeometry1D(UniformMeshGeometry):
    def __init__(self, origin, dx, numberOfFaces, numberOfCells, scale):
        self.origin = origin
        self.dx = dx
        self.numberOfFaces = numberOfFaces
        self.numberOfCells = numberOfCells
        self.scale = scale

        self._scaledGeometry = UniformMeshScaledGeometry1D(self)

    
    """Geometry properties"""
    @property
    def faceAreas(self):
        return numerix.ones(self.numberOfFaces,'d')

    @property
    def faceCenters(self):
        return numerix.arange(self.numberOfFaces)[numerix.NewAxis, ...] * self.dx + self.origin

    @property
    def faceNormals(self):
        faceNormals = numerix.ones((1, self.numberOfFaces), 'd')
        # The left-most face has neighboring cells None and the left-most cell.
        # We must reverse the normal to make fluxes work correctly.
        if self.numberOfFaces > 0:
            faceNormals[...,0] *= -1
        return faceNormals

    @property
    def orientedFaceNormals(self):
        return self.faceNormals

    @property
    def faceCellToCellNormals(self):
        return self.faceNormals
        
    @property
    def cellVolumes(self):
        return numerix.ones(self.numberOfCells, 'd') * self.dx

    @property
    def cellCenters(self):
        return ((numerix.arange(self.numberOfCells)[numerix.NewAxis, ...] + 0.5) \
                * self.dx + self.origin) * self.scale['length']

    @property
    def cellDistances(self):
        distances = numerix.ones(self.numberOfFaces, 'd')
        distances *= self.dx
        if len(distances) > 0:
            distances[0] = self.dx / 2.
            distances[-1] = self.dx / 2.
        return distances

    @property
    def faceTangents1(self):
        return numerix.zeros(self.numberOfFaces, 'd')[numerix.NewAxis, ...]

    @property
    def faceTangents2(self):
        return numerix.zeros(self.numberOfFaces, 'd')[numerix.NewAxis, ...]
    
    @property
    def cellToCellDistances(self):
        distances = MA.zeros((2, self.numberOfCells), 'd')
        distances[:] = self.dx
        if self.numberOfCells > 0:
            distances[0,0] = self.dx / 2.
            distances[1,-1] = self.dx / 2.
        return distances

    @property
    def cellNormals(self):
        normals = numerix.ones((1, 2, self.numberOfCells), 'd')
        if self.numberOfCells > 0:
            normals[:,0] = -1
        return normals
        
    @property
    def cellAreas(self):
        return numerix.ones((2, self.numberOfCells), 'd')

    @property
    def cellAreaProjections(self):
        return MA.array(self._getCellNormals())

    @property
    def faceCenters(self):
        return numerix.arange(self.numberOfFaces)[numerix.NewAxis, ...] * self.dx + self.origin
     
    @property
    def faceCellToCellNormals(self):
        return self.faceNormals

    @property
    def cellVolumes(self):
        return numerix.ones(self.numberOfCells, 'd') * self.dx
     
class UniformMeshGeometry2D(UniformMeshGeometry):

    def __init__(self, mesh):
        self.mesh = mesh

        self._scaledGeometry = UniformMeshScaledGeometry2D(
                                 self,
                                 self.mesh.numberOfFaces, 
                                 self.mesh.dx, 
                                 self.mesh.dy, 
                                 self.mesh.nx, 
                                 self.mesh.numberOfHorizontalFaces)

    @property
    def faceAreas(self):
        faceAreas = numerix.zeros(self.mesh.numberOfFaces, 'd')
        faceAreas[:self.mesh.numberOfHorizontalFaces] = self.mesh.dx
        faceAreas[self.mesh.numberOfHorizontalFaces:] = self.mesh.dy
        return faceAreas

    @property
    def faceNormals(self):
        normals = numerix.zeros((2, self.mesh.numberOfFaces), 'd')

        normals[1, :self.mesh.numberOfHorizontalFaces] = 1
        normals[1, :self.mesh.nx] = -1

        normals[0, self.mesh.numberOfHorizontalFaces:] = 1
        if self.mesh.numberOfVerticalColumns > 0:
            normals[0, self.mesh.numberOfHorizontalFaces::self.mesh.numberOfVerticalColumns] = -1

        return normals

    @property
    def faceCellToCellNormals(self):
        return self.faceNormals
        
    @property
    def cellVolumes(self):
        return numerix.ones(self.mesh.numberOfCells, 'd') * self.mesh.dx * self.mesh.dy

    @property
    def cellCenters(self):
        centers = numerix.zeros((2, self.mesh.nx, self.mesh.ny), 'd')
        indices = numerix.indices((self.mesh.nx, self.mesh.ny))
        centers[0] = (indices[0] + 0.5) * self.mesh.dx
        centers[1] = (indices[1] + 0.5) * self.mesh.dy
        return centers.reshape((2, self.mesh.numberOfCells), 
                               order="FORTRAN") + self.mesh.origin

    @property
    def cellDistances(self):
        Hdis = numerix.repeat((self.mesh.dy,), self.mesh.numberOfHorizontalFaces)
        Hdis = numerix.reshape(Hdis, (self.mesh.nx, self.mesh.numberOfHorizontalRows))
        if self.mesh.numberOfHorizontalRows > 0:
            Hdis[...,0] = self.mesh.dy / 2.
            Hdis[...,-1] = self.mesh.dy / 2.
        
        Vdis = numerix.repeat((self.mesh.dx,), self.mesh.numberOfFaces - self.mesh.numberOfHorizontalFaces)
        Vdis = numerix.reshape(Vdis, (self.mesh.numberOfVerticalColumns, self.mesh.ny))
        if self.mesh.numberOfVerticalColumns > 0:
            Vdis[0,...] = self.mesh.dx / 2.
            Vdis[-1,...] = self.mesh.dx / 2.

        return numerix.concatenate((numerix.reshape(numerix.swapaxes(Hdis,0,1), (self.mesh.numberOfHorizontalFaces,)), 
                                    numerix.reshape(numerix.swapaxes(Vdis,0,1), (self.mesh.numberOfFaces - self.mesh.numberOfHorizontalFaces,))))

    @property
    def faceToCellDistanceRatio(self):
        faceToCellDistanceRatios = numerix.zeros(self.mesh.numberOfFaces, 'd')
        faceToCellDistanceRatios[:] = 0.5
        faceToCellDistanceRatios[:self.mesh.nx] = 1.
        faceToCellDistanceRatios[self.mesh.numberOfHorizontalFaces - self.mesh.nx:self.mesh.numberOfHorizontalFaces] = 1.
        if self.mesh.numberOfVerticalColumns > 0:
            faceToCellDistanceRatios[self.mesh.numberOfHorizontalFaces::self.mesh.numberOfVerticalColumns] = 1.
            faceToCellDistanceRatios[(self.mesh.numberOfHorizontalFaces + self.mesh.nx)::self.mesh.numberOfVerticalColumns] = 1.
        return faceToCellDistanceRatios

    @property
    def faceToCellDistances(self):
        faceToCellDistances = numerix.zeros((2, self.mesh.numberOfFaces), 'd')
        distances = self.cellDistances
        ratios = self.faceToCellDistanceRatio
        faceToCellDistances[0] = distances * ratios
        faceToCellDistances[1] = distances * (1 - ratios)
        return faceToCellDistances
    
    @property
    def faceTangents1(self):
        tangents = numerix.zeros((2,self.mesh.numberOfFaces), 'd')

        if self.mesh.numberOfFaces > 0:
            tangents[0, :self.mesh.numberOfHorizontalFaces] = -1
            tangents[0, :self.mesh.nx] = 1        
            tangents[1, self.mesh.numberOfHorizontalFaces:] = 1
            tangents[1, self.mesh.numberOfHorizontalFaces::self.mesh.numberOfVerticalColumns] = -1

        return tangents
        
    @property
    def faceTangents2(self):
        return numerix.zeros((2, self.mesh.numberOfFaces), 'd')
    
    @property
    def cellToCellDistances(self):
        distances = numerix.zeros((4, self.mesh.nx, self.mesh.ny), 'd')
        distances[0] = self.mesh.dy
        distances[1] = self.mesh.dx
        distances[2] = self.mesh.dy
        distances[3] = self.mesh.dx
        
        if self.mesh.ny > 0:
            distances[0,..., 0] = self.mesh.dy / 2.
            distances[2,...,-1] = self.mesh.dy / 2.
        if self.mesh.nx > 0:
            distances[3, 0,...] = self.mesh.dx / 2.
            distances[1,-1,...] = self.mesh.dx / 2.
        
        return distances.reshape((4, self.mesh.numberOfCells), order="FORTRAN")


    @property
    def cellNormals(self):
        normals = numerix.zeros((2, 4, self.mesh.numberOfCells), 'd')
        normals[:, 0] = [[ 0], [-1]]
        normals[:, 1] = [[ 1], [ 0]]
        normals[:, 2] = [[ 0], [ 1]]
        normals[:, 3] = [[-1], [ 0]]

        return normals
        
    @property
    def cellAreas(self):
        areas = numerix.ones((4, self.mesh.numberOfCells), 'd')
        areas[0] = self.mesh.dx
        areas[1] = self.mesh.dy
        areas[2] = self.mesh.dx
        areas[3] = self.mesh.dy
        return areas

    @property
    def cellAreaProjections(self):
        return self.cellAreas * self.cellNormals

    @property
    def faceCenters(self):
        Hcen = numerix.zeros((2, self.mesh.nx, self.mesh.numberOfHorizontalRows), 'd')
        indices = numerix.indices((self.mesh.nx, self.mesh.numberOfHorizontalRows))
        Hcen[0,...] = (indices[0] + 0.5) * self.mesh.dx
        Hcen[1,...] = indices[1] * self.mesh.dy
        
        Vcen = numerix.zeros((2, self.mesh.numberOfVerticalColumns, self.mesh.ny), 'd')
        indices = numerix.indices((self.mesh.numberOfVerticalColumns, self.mesh.ny))
        Vcen[0,...] = indices[0] * self.mesh.dx
        Vcen[1,...] = (indices[1] + 0.5) * self.mesh.dy
        
        return numerix.concatenate((Hcen.reshape((2, self.mesh.numberOfHorizontalFaces), order="FORTRAN"),
                                    Vcen.reshape((2,
                                        self.mesh.numberOfVerticalFaces),
                                        order="FORTRAN")), axis=1) + self.mesh.origin
     
