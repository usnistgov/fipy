#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - a finite volume PDE solver in Python
 #
 #  FILE: "vanLeerConvectionTerm.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #
 # ========================================================================
 # This software was developed by employees of the National Institute
 # of Standards and Technology, an agency of the Federal Government.
 # Pursuant to title 17 section 105 of the United States Code,
 # works of NIST employees are not subject to copyright
 # protection, and this software is considered to be in the public domain.
 # FiPy is an experimental system.  NIST assumes no responsibility whatsoever
 # for its use by other parties, and makes no guarantees, expressed
 # or implied, about its quality, reliability, or any other characteristic.
 # We would appreciate acknowledgement if the document is used.
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

"""
"""

__docformat__ = 'restructuredtext'

from fipy.terms.explicitUpwindConvectionTerm import ExplicitUpwindConvectionTerm
from fipy.tools import numerix

__all__ = ["VanLeerConvectionTerm"]

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
