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

from fipy.meshes.geometries.abstractGeometries import AbstractMeshGeometry
        
class AbstractUniformGridGeometry(AbstractMeshGeometry):

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

    """Geometry properties common to 1D, 2D, 3D"""
    @property
    def orientedFaceNormals(self):
        return self.faceNormals

    def _getFaceToCellDistances(self):
        return self._faceToCellDistances
         
    @property
    def faceCellToCellNormals(self):
        return self.faceNormals
         
    def _setFaceToCellDistances(self, v):
        self._faceToCellDistances = v
        self._scaledGeometry._setScaledValues()

    faceToCellDistances = property(_getFaceToCellDistances,
                                   _setFaceToCellDistances)
     
 
