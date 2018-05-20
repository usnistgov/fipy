#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "modCellToFaceVariable.py"
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

__all__ = []

from fipy.tools import numerix
from fipy.tools import inline
from fipy.variables.arithmeticCellToFaceVariable import _ArithmeticCellToFaceVariable

class _ModCellToFaceVariable(_ArithmeticCellToFaceVariable):
    def __init__(self, var, modIn):
        _ArithmeticCellToFaceVariable.__init__(self,var)
        self.modIn = modIn

    if inline.doInline:
        def  _calcValue_(self, alpha, id1, id2):
            val = self._array.copy()

            inline._runInline(self.modIn + """
            int ID1 = id1[i];
            int ID2 = id2[i];
            double cell2 = var[ID2];
            val[i] = mod(cell2 - var[ID1]) * alpha[i] + var[ID1];
            """,var = self.var.numericValue,
                val = val,
                alpha = alpha,
                id1 = id1, id2 = id2,
                ni = self.mesh.numberOfFaces)

            return self._makeValue(value = val)
