from __future__ import unicode_literals
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
        Multiply a term

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

        Positive of a `Term`.

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
