#!/usr/bin/env python

##
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "surfactantConvectionVariable.py"
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

__docformat__ = 'restructuredtext'

__all__ = ['SurfactantConvectionVariable']

from fipy.tools.numerix import MA
from fipy.tools import numerix

from fipy.tools import vector

from fipy.variables.faceVariable import FaceVariable

class SurfactantConvectionVariable(FaceVariable):
    """

    Convection coefficient for the `ConservativeSurfactantEquation`.
    The coeff only has a value for a negative `distanceVar`.

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
           >>> answer[0,7] = -1
           >>> print numerix.allclose(SurfactantConvectionVariable(distanceVar).globalValue, answer)
           True

        Change the dimensions:

           >>> mesh = Grid2D(nx = 3, ny = 1, dx = .5, dy = .25)
           >>> distanceVar = DistanceVariable(mesh, value = (-.25, .25, .75))
           >>> answer[0,7] = -.5
           >>> print numerix.allclose(SurfactantConvectionVariable(distanceVar).globalValue, answer)
           True

        Two dimensional example:

           >>> mesh = Grid2D(nx = 2, ny = 2, dx = 1., dy = 1.)
           >>> distanceVar = DistanceVariable(mesh, value = (-1.5, -.5, -.5, .5))
            >>> answer = FaceVariable(mesh=mesh, rank=1, value=0.).globalValue
           >>> answer[1,2] = -.5
           >>> answer[1,3] = -1
           >>> answer[0,7] = -.5
           >>> answer[0,10] = -1
           >>> print numerix.allclose(SurfactantConvectionVariable(distanceVar).globalValue, answer)
           True

        Larger grid:

           >>> mesh = Grid2D(nx = 3, ny = 3, dx = 1., dy = 1.)
           >>> distanceVar = DistanceVariable(mesh, value = (1.5, .5 , 1.5,
           ...                                           .5 , -.5, .5 ,
           ...                                           1.5, .5 , 1.5))
            >>> answer = FaceVariable(mesh=mesh, rank=1, value=0.).globalValue
           >>> answer[1,4] = .25
           >>> answer[1,7] = -.25
           >>> answer[0,17] = .25
           >>> answer[0,18] = -.25
           >>> print numerix.allclose(SurfactantConvectionVariable(distanceVar).globalValue, answer)
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

        value = numerix.zeros((dim, Nfaces),'d')

        vector._putAdd(value, cellFaceIDs, alpha, mask=MA.getmask(MA.array(cellFaceIDs)))

##         value = numerix.reshape(value, (dim, Nfaces, dim))

        return -value / self.mesh._faceAreas

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
