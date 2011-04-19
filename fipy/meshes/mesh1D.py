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

"""Generic mesh class using numerix to do the calculations

    Meshes contain cells, faces, and vertices.

    This is built for a non-mixed element mesh.
"""

from fipy.tools import numerix
from fipy.tools.numerix import MA

from fipy.meshes.mesh import Mesh

from fipy.variables.faceVariable import FaceVariable

class Mesh1D(Mesh):
    
    def _calcScaleArea(self):
        return 1.

    def _calcScaleVolume(self):
        return self.scale['length']

    def _calcFaceAreas(self):
        return FaceVariable(mesh=self, value=1.)

    def _calcFaceNormals(self):
        faceNormals = FaceVariable(mesh=self, value=1., rank=1)
        # The left-most face has neighboring cells None and the left-most cell.
        # We must reverse the normal to make fluxes work correctly.
        if self.numberOfFaces > 0:
            faceNormals[...,0] *= -1
        return faceNormals

    def _calcFaceTangents(self):
        faceTangents1 = FaceVariable(mesh=self, value=0., rank=1)
        faceTangents2 = FaceVariable(mesh=self, value=0., rank=1)
        return faceTangents1, faceTangents2    

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
        
    @property
    def _VTKCellType(self):
        from enthought.tvtk.api import tvtk
        return tvtk.Line().cell_type

