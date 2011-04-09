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
from fipy.tools.numerix import MA

from fipy.variables.faceVariable import FaceVariable

from meshGeometry import MeshGeometry
from meshGeometry import ScaledMeshGeometry

class ScaledMeshGeometry1D(ScaledMeshGeometry):

    def _calcScaleArea(self):
        return 1.

    def _calcScaleVolume(self):
        return self.scale['length']
     
class MeshGeometry1D(MeshGeometry):
     
    def __init__(self, mesh, numberOfFaces, *args, **kwargs):
        self.numberOfFaces = numberOfFaces

        kwargs['ScaledGeom'] = ScaledMeshGeometry1D
        super(MeshGeometry1D, self).__init__(mesh, *args, **kwargs)

    def _calcFaceAreas(self):
        return FaceVariable(mesh=self.mesh, value=1.)

    def _calcFaceNormals(self):
        faceNormals = FaceVariable(mesh=self.mesh, value=1., rank=1)
        # The left-most face has neighboring cells None and the left-most cell.
        # We must reverse the normal to make fluxes work correctly.
        if self.numberOfFaces > 0:
            faceNormals[...,0] *= -1
        return faceNormals

    def _calcFaceTangents(self):
        faceTangents1 = FaceVariable(mesh=self.mesh, value=0., rank=1)
        faceTangents2 = FaceVariable(mesh=self.mesh, value=0., rank=1)
        return faceTangents1, faceTangents2
     
if __name__ == "__main__":
    import doctest
    doctest.testmod()

