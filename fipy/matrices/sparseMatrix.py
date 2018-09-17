#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "sparseMatrix.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #  Author: Maxsim Gibiansky <maxsim.gibiansky@nist.gov>
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
 #  2006-06-12 MLG 1.0 made abstract
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

__all__ = []

from fipy.tools import numerix

class _SparseMatrix(object):

    """
    .. attention:: This class is abstract. Always create one of its subclasses.
    """

    def __init__(self, mesh=None, bandwidth=0, matrix=None, sizeHint=None):
        pass

    matrix     = None
    numpyArray = property()
    _shape     = property()

    __array_priority__ = 100.0

    def __array_wrap(self, arr, context=None):
        if context is None:
            return arr
        else:
            return NotImplemented

    def copy(self):
        pass

    def __getitem__(self, index):
        pass

    def __str__(self):
        s = ''
        cellWidth = 11
        Irange, Jrange = self._range

        for j in Jrange:
            for i in Irange:
                v = self[j,i]
                if v == 0:
                    s += "---".center(cellWidth)
                else:
                    exp = numerix.log(abs(v))
                    if abs(exp) <= 4:
                        if exp < 0:
                            s += ("%9.6f" % v).ljust(cellWidth)
                        else:
                            s += ("%9.*f" % (6,v)).ljust(cellWidth)
                    else:
                        s += ("%9.2e" % v).ljust(cellWidth)
            s += "\n"
        return s[:-1]

    def __repr__(self):
        return repr(self.matrix)

    def __setitem__(self, index, value):
        pass

    def __add__(self, other):
        pass

    __radd__ = __add__

    def __iadd__(self, other):
        pass

    def __sub__(self, other):
        pass

    # Ask about this rsub
    def __rsub__(self, other):
        return -(__sub__(self, other))


    def __isub__(self, other):
        pass

    def __mul__(self, other):
        pass

    def __rmul__(self, other):
        pass

    def __neg__(self):
        return self * -1

    def __pos__(self):
        return self

##     def __eq__(self,other):
## 	return self.matrix.__eq__(other.matrix)

##     def transpose(self):
##         pass

    def put(self, vector, id1, id2):
        pass

    def putDiagonal(self, vector):
        pass

    def take(self, id1, id2):
        pass

    def takeDiagonal(self):
        pass

    def addAt(self, vector, id1, id2):
        pass

    def addAtDiagonal(self, vector):
        pass

    def exportMmf(self, filename):
        pass

##     def __array__(self):
##      shape = self._shape
##      indices = numerix.indices(shape)
##         numMatrix = self.take(indices[0].ravel(), indices[1].ravel())
##      return numerix.reshape(numMatrix, shape)
