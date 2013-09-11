#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "uniformGrid.py"
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
 #  
 # ###################################################################
 ##
__docformat__ = 'restructuredtext'
 
from fipy.meshes.abstractMesh import AbstractMesh

__all__ = ["UniformGrid"]

class UniformGrid(AbstractMesh):
    """Wrapped scaled geometry properties"""
    @property
    def _scaledFaceAreas(self):
        return self._faceAreas

    @property
    def _scaledCellVolumes(self):
        return self._cellVolumes

    @property
    def _scaledCellCenters(self):
        return self._cellCenters

    @property
    def _scaledCellDistances(self):
        return self._cellDistances

    @property
    def _scaledCellToCellDistances(self):
        return self._cellToCellDistances

    @property
    def _scaledFaceToCellDistances(self):
        return self._faceToCellDistances

    """Geometry properties common to 1D, 2D, 3D"""
    @property
    def _orientedFaceNormals(self):
        return self.faceNormals

    @property
    def _faceCellToCellNormals(self):
        return self.faceNormals
 
    def _getFaceToCellDistances(self):
        return self._internalFaceToCellDistances
                      
    def _setFaceToCellDistances(self, v):
        self._internalFaceToCellDistances = v
        self._setScaledValues()

    _faceToCellDistances = property(_getFaceToCellDistances,
                                    _setFaceToCellDistances)
