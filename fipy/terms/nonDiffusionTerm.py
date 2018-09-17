#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "nonDiffusionTerm.py"
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
from fipy.terms.unaryTerm import _UnaryTerm
from fipy.terms import TermMultiplyError

class _NonDiffusionTerm(_UnaryTerm):

    def __neg__(self):
        r"""
         Negate a `Term`.

           >>> -__NonDiffusionTerm(coeff=1.)
           __NonDiffusionTerm(coeff=-1.0)

        """
        if isinstance(self.coeff, (tuple, list)):
            return self.__class__(coeff=-numerix.array(self.coeff), var=self.var)
        else:
            return self.__class__(coeff=-self.coeff, var=self.var)

    def __mul__(self, other):
        r"""
        Mutiply a term

            >>> 2. * __NonDiffusionTerm(coeff=0.5)
            __NonDiffusionTerm(coeff=1.0)

        Test for ticket:291.

            >>> from fipy import PowerLawConvectionTerm
            >>> PowerLawConvectionTerm(coeff=[[1], [0]]) * 1.0
            PowerLawConvectionTerm(coeff=array([[ 1.],
                   [ 0.]]))

        """

        if isinstance(other, (int, float)):
            if isinstance(self.coeff, (list, tuple)):
                coeff = numerix.array(self.coeff)
            else:
                coeff = self.coeff
            return self.__class__(coeff=other * coeff, var=self.var)
        else:
            raise TermMultiplyError

    __rmul__ = __mul__

    @property
    def _diffusionVars(self):
        return []

    def _getDiagonalSign(self, transientGeomCoeff=None, diffusionGeomCoeff=None):
        if transientGeomCoeff is not None and diffusionGeomCoeff is not None:
            diagonalSign = numerix.where(numerix.array(numerix.all(transientGeomCoeff == 0, axis=-1)),
                                         numerix.array(2 * numerix.all(diffusionGeomCoeff[0] <= 0, axis=-1) - 1),
                                         numerix.array(2 * numerix.all(transientGeomCoeff >= 0, axis=-1) - 1))
        elif transientGeomCoeff is not None:
            diagonalSign = 2 * numerix.all(transientGeomCoeff >= 0, axis=-1) - 1
        elif diffusionGeomCoeff is not None:
            diagonalSign = 2 * numerix.all(diffusionGeomCoeff[0] <= 0, axis=-1) - 1
        else:
            diagonalSign = 1

        return diagonalSign

    def _test(self):
        r"""
        Test stuff.

         Subtract a `Term` from a `Term`, number or variable.

           >>> __NonDiffusionTerm(coeff=1.) - 10.
           (__NonDiffusionTerm(coeff=1.0) + -10.0)
           >>> __NonDiffusionTerm(coeff=1.) - __NonDiffusionTerm(coeff=2.)
           (__NonDiffusionTerm(coeff=1.0) + __NonDiffusionTerm(coeff=-2.0))

         Subtract a `Term`, number or variable from a `Term`.

           >>> 10. - __NonDiffusionTerm(coeff=1.)
           (__NonDiffusionTerm(coeff=-1.0) + 10.0)

        Add a `Term` to another `Term`, number or variable.

           >>> __NonDiffusionTerm(coeff=1.) + 10.
           (__NonDiffusionTerm(coeff=1.0) + 10.0)
           >>> __NonDiffusionTerm(coeff=1.) + __NonDiffusionTerm(coeff=2.)
           (__NonDiffusionTerm(coeff=1.0) + __NonDiffusionTerm(coeff=2.0))
           >>> 10. + __NonDiffusionTerm(coeff=1.)
           (__NonDiffusionTerm(coeff=1.0) + 10.0)

        Posate a `Term`.

           >>> +__NonDiffusionTerm(coeff=1.)
           __NonDiffusionTerm(coeff=1.0)

        This method allows `Terms` to be equated in a natural way. Note that the
        following does not return `False.`

           >>> __NonDiffusionTerm(coeff=1.) == __NonDiffusionTerm(coeff=2.)
           (__NonDiffusionTerm(coeff=1.0) + __NonDiffusionTerm(coeff=-2.0))

        it is equivalent to,

           >>> __NonDiffusionTerm(coeff=1.) - __NonDiffusionTerm(coeff=2.)
           (__NonDiffusionTerm(coeff=1.0) + __NonDiffusionTerm(coeff=-2.0))

        A `Term` can also equate with a number.

           >>> __NonDiffusionTerm(coeff=1.) == 1.
           (__NonDiffusionTerm(coeff=1.0) + -1.0)

        Likewise for integers.

           >>> __NonDiffusionTerm(coeff=1.) == 1
           (__NonDiffusionTerm(coeff=1.0) + -1)

        Equating to zero is allowed, of course

            >>> __NonDiffusionTerm(coeff=1.) == 0
            __NonDiffusionTerm(coeff=1.0)
            >>> 0 == __NonDiffusionTerm(coeff=1.)
            __NonDiffusionTerm(coeff=1.0)

        Divide a term

            >>> __NonDiffusionTerm(2.) / 2.
            __NonDiffusionTerm(coeff=1.0)

        Combine this equation with another

            >>> from fipy.variables.variable import Variable
            >>> eq1 = 10. + __NonDiffusionTerm(coeff=1., var=Variable(name='A'))
            >>> eq2 = 20. + __NonDiffusionTerm(coeff=2., var=Variable(name='B'))
            >>> eq1 & eq2
            ((__NonDiffusionTerm(coeff=1.0, var=A) + 10.0) & (__NonDiffusionTerm(coeff=2.0, var=B) + 20.0))

        """



class __NonDiffusionTerm(_NonDiffusionTerm):
    """
    Dummy subclass for tests
    """
    pass

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
