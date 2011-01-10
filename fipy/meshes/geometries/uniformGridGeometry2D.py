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
from fipy.tools import inline

from abstractGeometries import AbstractScaledMeshGeometry
from abstractUniformGeometries import AbstractUniformGridGeometry
 
class UniformGridScaledGeometry2D(AbstractScaledMeshGeometry):

    def __init__(self, geom, numFaces, dx, dy, nx, numHorizFaces):
        self.geom = geom
        self.numberOfFaces = numFaces
        self.dx = dx
        self.dy = dy
        self.nx = nx
        self.numberOfHorizontalFaces = numHorizFaces
      
    @property
    def orientedAreaProjections(self):
        return self.areaProjections

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

class UniformGridGeometry2D(AbstractUniformGridGeometry):

    def __init__(self, dx, dy,
                       nx, ny,
                       origin,
                       numberOfFaces,
                       numberOfHorizontalFaces,
                       numberOfVerticalFaces,
                       numberOfHorizontalRows,
                       numberOfVerticalColumns,
                       numberOfCells,
                       UniformScaledGeom=UniformGridScaledGeometry2D):

        self.dx = dx
        self.dy = dy
        self.nx = nx
        self.ny = ny
        self.origin = origin
        self.numberOfFaces = numberOfFaces
        self.numberOfHorizontalFaces = numberOfHorizontalFaces
        self.numberOfVerticalFaces = numberOfVerticalFaces
        self.numberOfHorizontalRows = numberOfHorizontalRows
        self.numberOfVerticalColumns = numberOfVerticalColumns
        self.numberOfCells = numberOfCells

        self._scaledGeometry = UniformScaledGeom(
                                 self,
                                 self.numberOfFaces, 
                                 self.dx, 
                                 self.dy, 
                                 self.nx, 
                                 self.numberOfHorizontalFaces)

    @property
    def faceAreas(self):
        faceAreas = numerix.zeros(self.numberOfFaces, 'd')
        faceAreas[:self.numberOfHorizontalFaces] = self.dx
        faceAreas[self.numberOfHorizontalFaces:] = self.dy
        return faceAreas

    @property
    def faceNormals(self):
        normals = numerix.zeros((2, self.numberOfFaces), 'd')

        normals[1, :self.numberOfHorizontalFaces] = 1
        normals[1, :self.nx] = -1

        normals[0, self.numberOfHorizontalFaces:] = 1
        if self.numberOfVerticalColumns > 0:
            normals[0, self.numberOfHorizontalFaces::self.numberOfVerticalColumns] = -1

        return normals

    @property
    def faceCellToCellNormals(self):
        return self.faceNormals
        
    @property
    def cellVolumes(self):
        return numerix.ones(self.numberOfCells, 'd') * self.dx * self.dy

    @property
    def cellCenters(self):
        centers = numerix.zeros((2, self.nx, self.ny), 'd')
        indices = numerix.indices((self.nx, self.ny))
        centers[0] = (indices[0] + 0.5) * self.dx
        centers[1] = (indices[1] + 0.5) * self.dy
        ccs = centers.reshape((2, self.numberOfCells), 
                               order="FORTRAN") + self.origin
        return ccs

    @property
    def cellDistances(self):
        Hdis = numerix.repeat((self.dy,), self.numberOfHorizontalFaces)
        Hdis = numerix.reshape(Hdis, (self.nx, self.numberOfHorizontalRows))
        if self.numberOfHorizontalRows > 0:
            Hdis[...,0] = self.dy / 2.
            Hdis[...,-1] = self.dy / 2.
        
        Vdis = numerix.repeat((self.dx,), self.numberOfFaces - self.numberOfHorizontalFaces)
        Vdis = numerix.reshape(Vdis, (self.numberOfVerticalColumns, self.ny))
        if self.numberOfVerticalColumns > 0:
            Vdis[0,...] = self.dx / 2.
            Vdis[-1,...] = self.dx / 2.

        return numerix.concatenate((numerix.reshape(numerix.swapaxes(Hdis,0,1), (self.numberOfHorizontalFaces,)), 
                                    numerix.reshape(numerix.swapaxes(Vdis,0,1), (self.numberOfFaces - self.numberOfHorizontalFaces,))))

    @property
    def faceToCellDistanceRatio(self):
        faceToCellDistanceRatios = numerix.zeros(self.numberOfFaces, 'd')
        faceToCellDistanceRatios[:] = 0.5
        faceToCellDistanceRatios[:self.nx] = 1.
        faceToCellDistanceRatios[self.numberOfHorizontalFaces - self.nx:self.numberOfHorizontalFaces] = 1.
        if self.numberOfVerticalColumns > 0:
            faceToCellDistanceRatios[self.numberOfHorizontalFaces::self.numberOfVerticalColumns] = 1.
            faceToCellDistanceRatios[(self.numberOfHorizontalFaces + self.nx)::self.numberOfVerticalColumns] = 1.
        return faceToCellDistanceRatios

    def _getFaceToCellDistances(self):
        if hasattr(self, "_faceToCellDistances"):
            """faces have been connected."""
            return self._faceToCellDistances

        else:
            faceToCellDistances = numerix.zeros((2, self.numberOfFaces), 'd')
            distances = self.cellDistances
            ratios = self.faceToCellDistanceRatio
            faceToCellDistances[0] = distances * ratios
            faceToCellDistances[1] = distances * (1 - ratios)
            return faceToCellDistances
    
    def _setFaceToCellDistances(self, v):
        """Exists only to allow `_connectFaces`."""
        self._faceToCellDistances = v

    faceToCellDistances = property(_getFaceToCellDistances,
                                   _setFaceToCellDistances)
     
    @property
    def faceTangents1(self):
        tangents = numerix.zeros((2,self.numberOfFaces), 'd')

        if self.numberOfFaces > 0:
            tangents[0, :self.numberOfHorizontalFaces] = -1
            tangents[0, :self.nx] = 1        
            tangents[1, self.numberOfHorizontalFaces:] = 1
            tangents[1, self.numberOfHorizontalFaces::self.numberOfVerticalColumns] = -1

        return tangents
        
    @property
    def faceTangents2(self):
        return numerix.zeros((2, self.numberOfFaces), 'd')
    
    @property
    def cellToCellDistances(self):
        distances = numerix.zeros((4, self.nx, self.ny), 'd')
        distances[0] = self.dy
        distances[1] = self.dx
        distances[2] = self.dy
        distances[3] = self.dx
        
        if self.ny > 0:
            distances[0,..., 0] = self.dy / 2.
            distances[2,...,-1] = self.dy / 2.
        if self.nx > 0:
            distances[3, 0,...] = self.dx / 2.
            distances[1,-1,...] = self.dx / 2.
        
        return distances.reshape((4, self.numberOfCells), order="FORTRAN")


    @property
    def cellNormals(self):
        normals = numerix.zeros((2, 4, self.numberOfCells), 'd')
        normals[:, 0] = [[ 0], [-1]]
        normals[:, 1] = [[ 1], [ 0]]
        normals[:, 2] = [[ 0], [ 1]]
        normals[:, 3] = [[-1], [ 0]]

        return normals
        
    @property
    def cellAreas(self):
        areas = numerix.ones((4, self.numberOfCells), 'd')
        areas[0] = self.dx
        areas[1] = self.dy
        areas[2] = self.dx
        areas[3] = self.dy
        return areas

    @property
    def cellAreaProjections(self):
        return self.cellAreas * self.cellNormals

    @property
    def faceCenters(self):
        Hcen = numerix.zeros((2, self.nx, self.numberOfHorizontalRows), 'd')
        indices = numerix.indices((self.nx, self.numberOfHorizontalRows))
        Hcen[0,...] = (indices[0] + 0.5) * self.dx
        Hcen[1,...] = indices[1] * self.dy
        
        Vcen = numerix.zeros((2, self.numberOfVerticalColumns, self.ny), 'd')
        indices = numerix.indices((self.numberOfVerticalColumns, self.ny))
        Vcen[0,...] = indices[0] * self.dx
        Vcen[1,...] = (indices[1] + 0.5) * self.dy
        
        return numerix.concatenate((Hcen.reshape((2, self.numberOfHorizontalFaces), order="FORTRAN"),
                                    Vcen.reshape((2,
                                        self.numberOfVerticalFaces),
                                        order="FORTRAN")), axis=1) + self.origin
  
