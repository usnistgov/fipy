#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "faceVariable.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
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

from fipy.variables.meshVariable import _MeshVariable
from fipy.tools import numerix

__all__ = ["FaceVariable"]

class FaceVariable(_MeshVariable):
    @property
    def _variableClass(self):
        return FaceVariable

    def _getShapeFromMesh(mesh):
        """
        Return the shape of this variable type, given a particular mesh.
        """
        return (mesh.numberOfFaces,)
    _getShapeFromMesh = staticmethod(_getShapeFromMesh)

    def _getArithmeticBaseClass(self, other = None):
        """
        Given `self` and `other`, return the desired base class for an operation
        result.
        """
        if other is None:
            return FaceVariable

        return _MeshVariable._getArithmeticBaseClass(self, other)

    def copy(self):
        return self._getArithmeticBaseClass()(mesh=self.mesh,
                                              name=self.name + "_copy",
                                              value=self.value)

    @property
    def globalValue(self):
        return self._getGlobalValue(self.mesh._localNonOverlappingFaceIDs,
                                    self.mesh._globalNonOverlappingFaceIDs)

    def setValue(self, value, unit = None, where = None):
        _MeshVariable.setValue(self, value=self._globalToLocalValue(value), unit=unit, where=where)

    @property
    def divergence(self):
        r"""the divergence of `self`, :math:`\vec{u}`,

        .. math:: \nabla\cdot\vec{u} \approx \frac{\sum_f (\vec{u}\cdot\hat{n})_f A_f}{V_P}

        Returns
        -------
        divergence : CellVariable
            one rank lower than `self`

        Examples
        --------

        >>> from fipy.meshes import Grid2D
        >>> from fipy.variables.cellVariable import CellVariable
        >>> mesh = Grid2D(nx=3, ny=2)
        >>> var = CellVariable(mesh=mesh, value=range(3*2))
        >>> print var.faceGrad.divergence
        [ 4.  3.  2. -2. -3. -4.]

        """
        if not hasattr(self, '_divergence'):
            from fipy.variables.addOverFacesVariable import _AddOverFacesVariable

            s = (slice(0,None,None),) + (numerix.newaxis,) * (len(self.shape) - 2) + (slice(0,None,None),)
            self._divergence = _AddOverFacesVariable((self * self.mesh._orientedAreaProjections[s]).sum(0))

        return self._divergence

    @property
    def _globalNumberOfElements(self):
        return self.mesh.globalNumberOfFaces

    @property
    def _globalOverlappingIDs(self):
        return self.mesh._globalOverlappingFaceIDs

    @property
    def _localNonOverlappingIDs(self):
        return self.mesh._localNonOverlappingFaceIDs

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
