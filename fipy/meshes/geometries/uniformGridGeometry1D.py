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
 
class UniformGridScaledGeometry1D(AbstractScaledMeshGeometry):

    def __init__(self, geom):
        self._geom = geom
    
    @property
    def faceToCellDistanceRatio(self):
        distances = numerix.ones(self._geom.numberOfFaces, 'd')
        distances *= 0.5
        if len(distances) > 0:
            distances[0] = 1
            distances[-1] = 1
        return distances

    @property
    def areaProjections(self):
        return self._geom.faceNormals
        
    @property
    def orientedAreaProjections(self):
        return self.areaProjections
        
    @property
    def _getFaceAspectRatios(self):
        return 1. / self._geom.cellDistances
     
class UniformGridGeometry1D(AbstractUniformGridGeometry):
    def __init__(self, mesh,
                       origin, 
                       dx, 
                       numberOfFaces, 
                       numberOfCells, 
                       scale,
                       ScaledGeomClass=UniformGridScaledGeometry1D):
        """TODO: Refactor args."""

        self.mesh = mesh
        self.origin = origin
        self.dx = dx
        self.numberOfFaces = numberOfFaces
        self.numberOfCells = numberOfCells
        self.scale = scale

        self._scaledGeometry = ScaledGeomClass(self)
    
    """Geometry properties"""
    @property
    def faceAreas(self):
        return FaceVariable(mesh=self.mesh, value=1.)

    @property
    def faceCenters(self):
        return numerix.arange(self.numberOfFaces)[numerix.NewAxis, ...] * self.dx + self.origin

    @property
    def faceNormals(self):
        faceNormals = FaceVariable(mesh=self.mesh, value=1., rank=1)
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
        return CellVariable(mesh=self.mesh, value=self.dx)

    @property
    def cellCenters(self):
        return CellVariable(mesh=self.mesh,
                            value=((numerix.arange(self.numberOfCells)[numerix.newaxis, ...] + 0.5) 
                                   * self.dx + self.origin) * self.scale['length'],
                            rank=1)

    @property
    def cellDistances(self):
        distances = FaceVariable(mesh=self.mesh, value=0.)
        distances[1:-1] = self.dx
        if len(distances) > 0:
            distances[0] = self.dx / 2.
            distances[-1] = self.dx / 2.
        return distances

    @property
    def faceTangents1(self):
        return FaceVariable(mesh=self.mesh, value=0., rank=1)

    @property
    def faceTangents2(self):
        return FaceVariable(mesh=self.mesh, value=0., rank=1)
    
    @property
    def cellToCellDistances(self):
        distances = MA.zeros((2, self.numberOfCells), 'd')
        distances[:] = self.dx
        if self.numberOfCells > 0:
            distances[0,0] = self.dx / 2.
            distances[1,-1] = self.dx / 2.
        return CellVariable(mesh=self.mesh, value=distances, elementshape=(2,))

    @property
    def cellNormals(self):
        normals = CellVariable(mesh=self.mesh, value=1., elementshape=(1,2))
        if self.numberOfCells > 0:
            normals[:,0] = -1
        return normals
        
    @property
    def cellAreas(self):
        return CellVariable(mesh=self.mesh, value=1, elementshape=(2,))

    @property
    def cellAreaProjections(self):
        return CellVariable(mesh=self.mesh,
                            value=MA.array(self.cellNormals.value),
                            elementshape=(1,2))

    @property
    def faceCenters(self):
        return FaceVariable(mesh=self.mesh,
                            value=numerix.arange(self.numberOfFaces)[numerix.newaxis, ...] * self.dx + self.origin,
                            rank=1)
     
    @property
    def faceCellToCellNormals(self):
        return self.faceNormals

    @property
    def cellVolumes(self):
        return numerix.ones(self.numberOfCells, 'd') * self.dx
    
