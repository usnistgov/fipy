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
 # of Standards and Technology, an agency of the Federal Government.
 # Pursuant to title 17 section 105 of the United States Code,
 # United States Code this software is not subject to copyright
 # protection, and this software is considered to be in the public domain.
 # FiPy is an experimental system.
 # NIST assumes no responsibility whatsoever for its use by whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 #
 # To the extent that NIST may hold copyright in countries other than the
 # United States, you are hereby granted the non-exclusive irrevocable and
 # unconditional right to print, publish, prepare derivative works and
 # distribute this software, in any medium, or authorize others to do so on
 # your behalf, on a royalty-free basis throughout the world.
 #
 # You may improve, modify, and create derivative works of the software or
 # any portion of the software, and you may copy and distribute such
 # modifications or works.  Modified works should carry a notice stating
 # that you changed the software and should note the date and nature of any
 # such change.  Please explicitly acknowledge the National Institute of
 # Standards and Technology as the original source.
 #
 # This software can be redistributed and/or modified freely provided that
 # any derivative works bear some notice that they are derived from it, and
 # any modified versions bear some notice that they have been modified.
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
