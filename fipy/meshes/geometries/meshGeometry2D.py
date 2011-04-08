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

from meshGeometry import MeshGeometry
from meshGeometry import ScaledMeshGeometry

class ScaledMeshGeometry2D(ScaledMeshGeometry):
    
    def _calcScaleArea(self):
        return self.scale['length']

    def _calcScaleVolume(self):
        return self.scale['length']**2
     
class MeshGeometry2D(MeshGeometry):
     
    def __init__(self, *args, **kwargs):
        kwargs['ScaledGeom'] = ScaledMeshGeometry2D
        super(MeshGeometry2D, self).__init__(*args, **kwargs)

    def _calcFaceAreas(self):
        faceVertexCoords = numerix.take(self.vertexCoords, self.faceVertexIDs, axis=1)
        tangent = faceVertexCoords[:,1] - faceVertexCoords[:,0]
        return numerix.sqrtDot(tangent, tangent)

    def _calcFaceNormals(self):
        faceVertexCoords = numerix.take(self.vertexCoords, self.faceVertexIDs, axis=1)
        t1 = faceVertexCoords[:,1,:] - faceVertexCoords[:,0,:]
        rot = numerix.eye(N=self.dim)
        rot[0:2,0:2] = numerix.array((( 0, 1),
                                      (-1, 0)))
        faceNormals = t1.dot(rot) / t1.mag
        
        orientation = 1 - 2 * (numerix.dot(faceNormals, self.cellDistanceVectors) < 0)
        return faceNormals * orientation

    def _calcFaceTangents(self):
        rot = numerix.eye(N=self.getDim())
        rot[0:2,0:2] = numerix.array((( 0, 1),
                                      (-1, 0)))
        tmp = self.faceNormals.dot(rot)
        faceTangents1 = tmp / tmp.mag
        faceTangents2 = FaceVariable(mesh=self, value=0., rank=1)
        return faceTangents1, faceTangents2
    
