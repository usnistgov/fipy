from __future__ import division
from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy.terms.explicitUpwindConvectionTerm import ExplicitUpwindConvectionTerm
from fipy.tools import numerix

__all__ = ["VanLeerConvectionTerm"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class VanLeerConvectionTerm(ExplicitUpwindConvectionTerm):

    def _getGradient(self, normalGradient, gradUpwind):
        gradUpUpwind = -gradUpwind + 2 * normalGradient

        avg = 0.5 * (abs(gradUpwind) + abs(gradUpUpwind))
        min3 = numerix.minimum(numerix.minimum(abs(2 * gradUpwind),
                                               abs(2 * gradUpUpwind)), avg)

        grad = numerix.where(gradUpwind * gradUpUpwind < 0.,
                             0.,
                             numerix.where(gradUpUpwind > 0.,
                                           min3,
                                           -min3))

        return grad

    def _getOldAdjacentValues(self, oldArray, id1, id2, dt):
        oldArray1, oldArray2 = ExplicitUpwindConvectionTerm._getOldAdjacentValues(self, oldArray, id1, id2, dt)

        mesh = oldArray.mesh

        interiorIDs = numerix.nonzero(mesh.interiorFaces)[0]
        interiorFaceAreas = numerix.take(mesh._faceAreas, interiorIDs)
        interiorFaceNormals = numerix.take(mesh._orientedFaceNormals, interiorIDs, axis=-1)

        # Courant-Friedrichs-Levy number
        interiorCFL = abs(numerix.take(self._getGeomCoeff(oldArray), interiorIDs)) * dt

        gradUpwind = (oldArray2 - oldArray1) / numerix.take(mesh._cellDistances, interiorIDs)

        vol1 = numerix.take(mesh.cellVolumes, id1)

        oldArray1 += 0.5 * self._getGradient(numerix.dot(numerix.take(oldArray.grad, id1, axis=-1), interiorFaceNormals), gradUpwind) \
            * (vol1 - interiorCFL) / interiorFaceAreas

        vol2 = numerix.take(mesh.cellVolumes, id2)

        oldArray2 += 0.5 * self._getGradient(numerix.dot(numerix.take(oldArray.grad, id2, axis=-1), -interiorFaceNormals), -gradUpwind) \
            * (vol2 - interiorCFL) / interiorFaceAreas

        return oldArray1, oldArray2

    def _test(self):
        """
        Test for ticket:441.

        >>> from fipy import *
        >>> m = Grid1D()
        >>> c = CellVariable(mesh=m)
        >>> e = VanLeerConvectionTerm(((1,),))
        >>> e.solve(c) # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
           ...
        TransientTermError: The equation requires a TransientTerm with explicit convection.

        """

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
