#!/usr/bin/env python

## 
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "mesh1D.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
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

"""Generic mesh class using numerix to do the calculations

    Meshes contain cells, faces, and vertices.

    This is built for a non-mixed element mesh.
"""

from fipy.tools import numerix
from fipy.tools.numerix import MA

from fipy.meshes.numMesh.mesh import Mesh

class Mesh1D(Mesh):
    def _calcFaceAreas(self):
        self.faceAreas = numerix.ones(self.numberOfFaces, 'd')

    def _calcFaceNormals(self):
        self.faceNormals = numerix.array((numerix.ones(self.numberOfFaces, 'd'),))
        # The left-most face has neighboring cells None and the left-most cell.
        # We must reverse the normal to make fluxes work correctly.
        if self.numberOfFaces > 0:
            self.faceNormals[...,0] = -self.faceNormals[...,0]

    def _calcFaceTangents(self):
        self.faceTangents1 = numerix.zeros(self.numberOfFaces, 'd')[numerix.NewAxis, ...]
        self.faceTangents2 = numerix.zeros(self.numberOfFaces, 'd')[numerix.NewAxis, ...]

    def _calcHigherOrderScalings(self):
        self.scale['area'] = 1.
        self.scale['volume'] = self.scale['length']

    def _translate(self, vector):
        newCoords = self.vertexCoords + vector
        newmesh = Mesh1D(newCoords, numerix.array(self.faceVertexIDs), numerix.array(self.cellFaceIDs))
        return newmesh

    def __mul__(self, factor):
        newCoords = self.vertexCoords * factor
        newmesh = Mesh1D(newCoords, numerix.array(self.faceVertexIDs), numerix.array(self.cellFaceIDs))
        return newmesh

    @property
    def _concatenatedClass(self):
        return Mesh1D

    def _isOrthogonal(self):
        return True
        
    def _getVTKCellType(self):
        from enthought.tvtk.api import tvtk
        return tvtk.Line().cell_type

