#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "baseUpwindConvectionTerm.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
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

__all__ = []

from fipy.terms.abstractConvectionTerm import _AbstractConvectionTerm
from fipy.variables.faceVariable import FaceVariable
from fipy.tools.dimensions.physicalField import PhysicalField
from fipy.tools import inline
from fipy.tools import numerix

class _UpwindConvectionTermAlpha(FaceVariable):
    def __init__(self, P):
        FaceVariable.__init__(self, mesh=P.mesh, elementshape=P.shape[:-1])
        self.P = self._requires(P)

    if inline.doInline:
        def _calcValue(self):
            P  = self.P.numericValue
            alpha = self._array.copy()

            inline._runInline("""
                alpha[i] = 0.5;

                if (P[i] > 0.) {
                    alpha[i] = 1.;
                } else {
                    alpha[i] = 0.;
                }
            """,
            alpha=alpha, P=P,
            ni = len(P.flat))

            return self._makeValue(value=alpha)
    else:
        def _calcValue(self):
            P  = self.P.numericValue
            alpha = numerix.where(P > 0., 1., 0.)
            return PhysicalField(value=alpha)

class _AbstractUpwindConvectionTerm(_AbstractConvectionTerm):
    def _alpha(self, P):
        return _UpwindConvectionTermAlpha(P)

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
