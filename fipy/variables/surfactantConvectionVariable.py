from __future__ import division
from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

__all__ = ['SurfactantConvectionVariable']
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

from fipy.tools.numerix import MA
from fipy.tools import numerix

from fipy.tools import vector

from fipy.variables.faceVariable import FaceVariable

class SurfactantConvectionVariable(FaceVariable):
    """

    Convection coefficient for the `ConservativeSurfactantEquation`.
    The coefficient only has a value for a negative `distanceVar`.

    """

    def __init__(self, distanceVar):
        """

        Simple one dimensional test:


           >>> from fipy.variables.cellVariable import CellVariable
           >>> from fipy.meshes import Grid2D
           >>> mesh = Grid2D(nx = 3, ny = 1, dx = 1., dy = 1.)
           >>> from fipy.variables.distanceVariable import DistanceVariable
           >>> distanceVar = DistanceVariable(mesh, value = (-.5, .5, 1.5))
           >>> ## answer = numerix.zeros((2, mesh.numberOfFaces),'d')
           >>> answer = FaceVariable(mesh=mesh, rank=1, value=0.).globalValue
           >>> answer[0, 7] = -1
           >>> print(numerix.allclose(SurfactantConvectionVariable(distanceVar).globalValue, answer))
           True

        Change the dimensions:

           >>> mesh = Grid2D(nx = 3, ny = 1, dx = .5, dy = .25)
           >>> distanceVar = DistanceVariable(mesh, value = (-.25, .25, .75))
           >>> answer[0, 7] = -.5
           >>> print(numerix.allclose(SurfactantConvectionVariable(distanceVar).globalValue, answer))
           True

        Two dimensional example:

           >>> mesh = Grid2D(nx = 2, ny = 2, dx = 1., dy = 1.)
           >>> distanceVar = DistanceVariable(mesh, value = (-1.5, -.5, -.5, .5))
            >>> answer = FaceVariable(mesh=mesh, rank=1, value=0.).globalValue
           >>> answer[1, 2] = -.5
           >>> answer[1, 3] = -1
           >>> answer[0, 7] = -.5
           >>> answer[0, 10] = -1
           >>> print(numerix.allclose(SurfactantConvectionVariable(distanceVar).globalValue, answer))
           True

        Larger grid:

           >>> mesh = Grid2D(nx = 3, ny = 3, dx = 1., dy = 1.)
           >>> distanceVar = DistanceVariable(mesh, value = (1.5, .5, 1.5,
           ...                                           .5, -.5, .5,
           ...                                           1.5, .5, 1.5))
            >>> answer = FaceVariable(mesh=mesh, rank=1, value=0.).globalValue
           >>> answer[1, 4] = .25
           >>> answer[1, 7] = -.25
           >>> answer[0, 17] = .25
           >>> answer[0, 18] = -.25
           >>> print(numerix.allclose(SurfactantConvectionVariable(distanceVar).globalValue, answer))
           True

        """

        FaceVariable.__init__(self, mesh=distanceVar.mesh, name='surfactant convection', rank=1)
        self.distanceVar = self._requires(distanceVar)

    def _calcValue(self):

        Nfaces = self.mesh.numberOfFaces
        M = self.mesh._maxFacesPerCell
        dim = self.mesh.dim
        cellFaceIDs = self.mesh.cellFaceIDs

        faceNormalAreas = self.distanceVar._levelSetNormals * self.mesh._faceAreas

        cellFaceNormalAreas = numerix.array(MA.filled(numerix.take(faceNormalAreas, cellFaceIDs, axis=-1), 0))
        norms = numerix.array(MA.filled(MA.array(self.mesh._cellNormals), 0))

        alpha = numerix.dot(cellFaceNormalAreas, norms)
        alpha = numerix.where(alpha > 0, alpha, 0)

        alphasum = numerix.sum(alpha, axis=0)
        alphasum += (alphasum < 1e-100) * 1.0
        alpha = alpha / alphasum

        phi = numerix.repeat(self.distanceVar[numerix.newaxis, ...], M, axis=0)
        alpha = numerix.where(phi > 0., 0, alpha)

        volumes = numerix.array(self.mesh.cellVolumes)
        alpha = alpha * volumes * norms

        value = numerix.zeros((dim, Nfaces), 'd')

        vector._putAdd(value, cellFaceIDs, alpha, mask=MA.getmask(MA.array(cellFaceIDs)))

##         value = numerix.reshape(value, (dim, Nfaces, dim))

        return -value / self.mesh._faceAreas

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()


