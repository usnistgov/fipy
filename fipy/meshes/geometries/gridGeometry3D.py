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

from meshGeometry import MeshGeometry
from meshGeometry import ScaledMeshGeometry

class ScaledGridGeometry3D(ScaledMeshGeometry):
    
    def _calcScaleArea(self):
        return self.scale['length']**2

    def _calcScaleVolume(self):
        return self.scale['length']**3

class GridGeometry3D(MeshGeometry):

    def __init__(self, nx, ny, nz,
                       numberOfFaces,
                       numberOfXYFaces,
                       numberOfXZFaces,
                       numberOfYZFaces, 
                       *args,
                       **kwargs):

        self.nx = nx
        self.ny = ny
        self.nz = nz

        self.numberOfFaces = numberOfFaces
        self.numberOfXYFaces = numberOfXYFaces
        self.numberOfXZFaces = numberOfXZFaces
        self.numberOfYZFaces = numberOfYZFaces

        kwargs['ScaledGeom'] = ScaledGridGeometry3D
        super(GridGeometry3D, self).__init__(*args, **kwargs)

    def _calcFaceNormals(self):
        XYFaceNormals = numerix.zeros((3, self.numberOfXYFaces))
        XYFaceNormals[2, (self.nx * self.ny):] = 1
        XYFaceNormals[2, :(self.nx * self.ny)] = -1
        XZFaceNormals = numerix.zeros((3, self.numberOfXZFaces))
        xzd = numerix.arange(self.numberOfXZFaces)
        xzd = xzd % (self.nx * (self.ny + 1))
        xzd = (xzd < self.nx)
        xzd = 1 - (2 * xzd)
        XZFaceNormals[1, :] = xzd
        YZFaceNormals = numerix.zeros((3, self.numberOfYZFaces))
        YZFaceNormals[0, :] = 1
        YZFaceNormals[0, ::self.nx + 1] = -1
        return numerix.concatenate((XYFaceNormals, 
                                    XZFaceNormals, 
                                    YZFaceNormals), 
                                   axis=-1)
        
    def _calcFaceTangents(self):
        ## need to see whether order matters.
        faceTangents1 = numerix.zeros((3, self.numberOfFaces), 'd')
        faceTangents2 = numerix.zeros((3, self.numberOfFaces), 'd')
        ## XY faces
        faceTangents1[0, :self.numberOfXYFaces] = 1.
        faceTangents2[1, :self.numberOfXYFaces] = 1.
        ## XZ faces
        faceTangents1[0, self.numberOfXYFaces:self.numberOfXYFaces + self.numberOfXZFaces] = 1.
        faceTangents2[2, self.numberOfXYFaces:self.numberOfXYFaces + self.numberOfXZFaces] = 1.
        ## YZ faces
        faceTangents1[1, self.numberOfXYFaces + self.numberOfXZFaces:] = 1.
        faceTangents2[2, self.numberOfXYFaces + self.numberOfXZFaces:] = 1.
        return faceTangents1, faceTangents2
 
