#!/usr/bin/env python

## -*-Pyth-*-
 # #############################################################################
 # FiPy - a finite volume PDE solver in Python
 #
 # FILE: "binaryOperatorVariable.py"
 #
 # Author: Jonathan Guyer <guyer@nist.gov>
 # Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 # Author: James Warren   <jwarren@nist.gov>
 #   mail: NIST
 #    www: <http://www.ctcms.nist.gov/fipy/>
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
 # #############################################################################
 ##

__docformat__ = 'restructuredtext'

__all__ = []

from fipy.tools import numerix

def _BinaryOperatorVariable(operatorClass=None):
    """
    Test BinOp pickling

        >>> from fipy import Grid1D, FaceVariable, CellVariable, dump, Variable
        >>> import os, sys
        >>> m = Grid1D()
        >>> vs = (CellVariable(mesh=m, value=2.), FaceVariable(mesh=m, value=3.), Variable(4))
        >>> tmp = []
        >>> for v in vs:
        ...     (f, n) = dump.write(v * v)
        ...     tmp += [dump.read(n)]
        ...     if f is not None:
        ...         os.close(f)
        ...         os.remove(n)
        >>> for v in tmp:
        ...     print v.__class__
        <class 'fipy.variables.cellVariable.CellVariable'>
        <class 'fipy.variables.faceVariable.FaceVariable'>
        <class 'fipy.variables.variable.Variable'>
        >>> print tmp[0].allclose(4.)
        True
        >>> print tmp[1].allclose(9.)
        True
        >>> print tmp[2].allclose(16)
        True

    """
    # declare a binary operator class with the desired base class
    class binOp(operatorClass):

        def _calcValue_(self):
            from fipy.variables.variable import Variable
            if isinstance(self.var[1], Variable):
                val1 = self.var[1].value
            else:
                if type(self.var[1]) is type(''):
                    self.var[1] = physicalField.PhysicalField(value=self.var[1])
                val1 = self.var[1]

            return self.op(self.var[0].value, val1)

        @property
        def unit(self):
            if self._unit is None:
                try:
                    var = self._varProxy
                    return self._extractUnit(self.op(var[0], var[1]))
                except:
                    return self._extractUnit(self._calcValue_())
            else:
                return self._unit

        def _getRepresentation(self, style="__repr__", argDict={}, id=id, freshen=False):
            self.id = id
            if (style == "__repr__") and hasattr(self, '_name') and len(self._name) > 0:
                return self._name
            else:
                return "(" + operatorClass._getRepresentation(self, style=style, argDict=argDict, id=id, freshen=freshen) + ")"

    return binOp

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
