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

from abstractGeometries import AbstractScaledMeshGeometry
from abstractUniformGeometries import AbstractUniformGridGeometry
 
class UniformGridScaledGeometry1D(AbstractScaledMeshGeometry):

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
        return self.areaProjections
        
    @property
    def _getFaceAspectRatios(self):
        return 1. / self.geom.cellDistances
     
class UniformGridGeometry1D(AbstractUniformGridGeometry):
    def __init__(self, origin, 
                       dx, 
                       numberOfFaces, 
                       numberOfCells, 
                       scale,
                       ScaledGeomClass=UniformGridScaledGeometry1D):
        """TODO: Refactor args."""

        self.origin = origin
        self.dx = dx
        self.numberOfFaces = numberOfFaces
        self.numberOfCells = numberOfCells
        self.scale = scale

        self._scaledGeometry = ScaledGeomClass(self)
    
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
        ccs = ((numerix.arange(self.numberOfCells)[numerix.NewAxis, ...] + 0.5) \
               * self.dx + self.origin) * self.scale['length']
        return ccs

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
        return MA.array(self.cellNormals)

    @property
    def faceCenters(self):
        return numerix.arange(self.numberOfFaces)[numerix.NewAxis, ...] * self.dx + self.origin
     
    @property
    def faceCellToCellNormals(self):
        return self.faceNormals

    @property
    def cellVolumes(self):
        return numerix.ones(self.numberOfCells, 'd') * self.dx
    
