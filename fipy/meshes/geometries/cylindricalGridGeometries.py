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

from meshGeometry1D import MeshGeometry1D
from meshGeometry2D import MeshGeometry2D
from uniformMeshGeometry import UniformMeshGeometry1D
from uniformMeshGeometry import UniformMeshGeometry2D
from uniformMeshGeometry import UniformMeshScaledGeometry1D
from uniformMeshGeometry import UniformMeshScaledGeometry2D

class CylindricalGridGeometry1D(MeshGeometry1D):
     
    def __init__(self, mesh, scaleLength):
        self.origin = mesh.origin

        super(CylindricalGridGeometry1D, self).__init__(mesh, 
                                                        scaleLength)
     
    def _calcFaceCenters(self):
        faceCenters = super(CylindricalGridGeometry1D, self)._calcFaceCenters()
        return faceCenters + self.origin
    
    def _calcFaceAreas(self):
        return self._calcFaceCenters()[0]

    def _calcCellVolumes(self):
        return super(CylindricalGridGeometry1D, self)._calcCellVolumes() / 2.
     
class CylindricalGridGeometry2D(MeshGeometry2D):
     
    def __init__(self, mesh, scaleLength):
        self.origin = mesh.origin
        self.mesh = mesh

        super(CylindricalGridGeometry2D, self).__init__(mesh, 
                                                        scaleLength)
     
    @property
    def faceAreas(self):
        return self._faceAreas * self.faceCenters[0]

    @property
    def cellVolumes(self):
        return super(CylindricalGridGeometry2D, self).cellVolumes * self.cellCenters[0]
     
    @property
    def cellCenters(self):
        return super(CylindricalGridGeometry2D, self).cellCenters + self.origin

    @property
    def faceCenters(self):
        return self._calcFaceCenters() + self.origin
     
class CylindricalUniformGridScaledGeometry1D(UniformMeshScaledGeometry1D):

    @property
    def faceAspectRatios(self):
        return self.geom.faceAreas / self.geom.cellDistances
    
    @property
    def areaProjections(self):
        return self.geom.faceNormals * self.geom.faceAreas
       
class CylindricalUniformGridGeometry1D(UniformMeshGeometry1D):
    def __init__(self, mesh, origin, dx, numberOfFaces, numberOfCells, scale,
                 ScaledGeomClass=UniformMeshScaledGeometry1D):
        super(CylindricalUniformGridGeometry1D, self).__init__(
                mesh,
                origin,
                dx,
                numberOfFaces,
                numberOfCells,
                scale,
                ScaledGeomClass=CylindricalUniformGridScaledGeometry1D)
        
    @property
    def faceAreas(self):
        return self.faceCenters[0]

    @property
    def cellVolumes(self):
        return self.cellCenters[0] * self.dx

    @property
    def cellAreas(self):
        return numerix.array((self.faceAreas[:-1], self.faceAreas[1:]))

    @property
    def cellAreaProjections(self):
        return MA.array(self.cellNormals) * self.cellAreas

      
class CylindricalUniformGridScaledGeometry2D(UniformMeshScaledGeometry2D):

    def _calcAreaProjections(self):
        return self._getAreaProjectionsPy()
     
class CylindricalUniformGridGeometry2D(UniformMeshGeometry2D):

    def __init__(self, mesh):
        super(CylindricalUniformGridGeometry2D, self).__init__(mesh,
            UniformScaledGeom=CylindricalUniformGridScaledGeometry2D)

    @property
    def faceAreas(self):
        faceAreas = numerix.zeros(self.mesh.numberOfFaces, 'd')
        faceAreas[:self.mesh.numberOfHorizontalFaces] = self.mesh.dx
        faceAreas[self.mesh.numberOfHorizontalFaces:] = self.mesh.dy
        return faceAreas * self.faceCenters[0]
        
    @property
    def cellVolumes(self):
        return numerix.ones(self.mesh.numberOfCells, 'd') * self.mesh.dx \
                * self.mesh.dy * self.cellCenters[0]

    @property
    def cellAreas(self):
        areas = numerix.ones((4, self.mesh.numberOfCells), 'd')
        areas[0] = self.mesh.dx * self.cellCenters[0]
        areas[1] = self.mesh.dy * (self.cellCenters[0] + self.mesh.dx / 2)
        areas[2] = self.mesh.dx * self.cellCenters[0]
        areas[3] = self.mesh.dy * (self.cellCenters[0] - self.mesh.dx / 2)
        return areas
     
