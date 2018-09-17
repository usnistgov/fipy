#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "implicitSourceTerm.py"
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

from fipy.terms.sourceTerm import SourceTerm
from fipy.tools import numerix

__all__ = ["ImplicitSourceTerm"]

class ImplicitSourceTerm(SourceTerm):
    r"""

    The `ImplicitSourceTerm` represents

    .. math::

       \int_V \phi S \,dV \simeq \phi_P S_P V_P

    where :math:`S` is the `coeff` value.
    """

    def _getWeight(self, var, transientGeomCoeff=None, diffusionGeomCoeff=None):
        """
        Test for a bug due to the sign operator not being updating
        correctly.

            >>> from fipy import *
            >>> m = Grid1D(nx=1)
            >>> v = CellVariable(mesh=m, value=1.)
            >>> eq = TransientTerm() == ImplicitSourceTerm(v)
            >>> eq.solve(v, dt=1.)
            >>> print v
            [ 2.]
            >>> v.setValue(-1.)
            >>> eq.solve(v, dt=1.)
            >>> print v
            [-0.5]

        """

        coeff = self._getGeomCoeff(var)
        diagonalSign = self._getDiagonalSign(transientGeomCoeff, diffusionGeomCoeff)
        combinedSign = numerix.array(diagonalSign)[...,numerix.newaxis] * numerix.sign(coeff)

        return {'diagonal' : (combinedSign >= 0),
                'old value' : numerix.zeros(var.shape, 'd'),
                'b vector' :  -var * (combinedSign < 0),
                'new value' : numerix.zeros(var.shape, 'd')}

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
