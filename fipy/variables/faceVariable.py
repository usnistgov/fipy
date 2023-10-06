from __future__ import unicode_literals
from fipy.variables.meshVariable import MeshVariable
from fipy.tools import numerix

__all__ = ["FaceVariable"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class FaceVariable(MeshVariable):
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

        return MeshVariable._getArithmeticBaseClass(self, other)

    def copy(self):
        return self._getArithmeticBaseClass()(mesh=self.mesh,
                                              name=self.name + "_copy",
                                              value=self.value)

    @property
    def globalValue(self):
        ownedFaceIDs = self._localNonOverlappingIDs
        return self._getGlobalValue(ownedFaceIDs,
                                    self._globalOverlappingIDs[..., ownedFaceIDs])

    def setValue(self, value, unit = None, where = None):
        MeshVariable.setValue(self, value=self._globalToLocalValue(value), unit=unit, where=where)

    @property
    def divergence(self):
        r"""the divergence of `self`, :math:`\vec{u}`,

        .. math:: \nabla\cdot\vec{u} \approx \frac{\sum_f (\vec{u}\cdot\hat{n})_f A_f}{V_P}

        Returns
        -------
        divergence : fipy.variables.cellVariable.CellVariable
            one rank lower than `self`

        Examples
        --------

        >>> from fipy.meshes import Grid2D
        >>> from fipy.variables.cellVariable import CellVariable
        >>> mesh = Grid2D(nx=3, ny=2)
        >>> from builtins import range
        >>> var = CellVariable(mesh=mesh, value=list(range(3*2)))
        >>> print(var.faceGrad.divergence)
        [ 4.  3.  2. -2. -3. -4.]

        """
        if not hasattr(self, '_divergence'):
            from fipy.variables.addOverFacesVariable import _AddOverFacesVariable

            s = (slice(0, None, None),) + (numerix.newaxis,) * (len(self.shape) - 2) + (slice(0, None, None),)
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
        return self.mesh.topology._ownedFaceIDs

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
