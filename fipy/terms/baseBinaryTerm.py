#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - a finite volume PDE solver in Python
 # 
 #  FILE: "baseBinaryTerm.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #  
 # ========================================================================
 # This document was prepared at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this document is not subject to copyright
 # protection and is in the public domain.  summationTerm.py
 # is an experimental work.  NIST assumes no responsibility whatsoever
 # for its use by other parties, and makes no guarantees, expressed
 # or implied, about its quality, reliability, or any other characteristic.
 # We would appreciate acknowledgement if the document is used.
 # 
 # This document can be redistributed and/or modified freely
 # provided that any derivative works bear some notice that they are
 # derived from it, and any modified versions bear some notice that
 # they have been modified.
 # ========================================================================
 #  See the file "license.terms" for information on usage and  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #  
 # ###################################################################
 ##

import os

from fipy.terms.term import Term
from fipy.terms.explicitSourceTerm import _ExplicitSourceTerm

class _BaseBinaryTerm(Term):
    def __init__(self, term, other):

        if not isinstance(other, Term):
            other = _ExplicitSourceTerm(coeff=other, var=term.var)

        self.term = term
        self.other = other

        if term.var is None:
            if other.var is None:
                pass
            else:
                raise Exception, 'Terms with explicit Variables cannot mix with Terms with implicit Variables'
        else:
            if other.var is None:
                raise Exception, 'Terms with explicit Variables cannot mix with Terms with implicit Variables'

	Term.__init__(self, var=self._getVars()[0])

    def _getVars(self):
        
        if not hasattr(self, '_vars'):
            seen = set()
            seq = self.term._getVars() + self.other._getVars()
            self._vars = [x for x in seq if x not in seen and not seen.add(x)]
        return self._vars

    def _addNone(self, arg0, arg1):
        if arg0 is None and arg1 is None:
            return None
        elif arg0 is None:
            return arg1
        elif arg1 is None:
            return arg0
        else:
            return arg0 + arg1

    def _getTransientGeomCoeff(self, mesh):
        return self._addNone(self.term._getTransientGeomCoeff(mesh), self.other._getTransientGeomCoeff(mesh))

    def _getDiffusionGeomCoeff(self, mesh):
        return self._addNone(self.term._getDiffusionGeomCoeff(mesh), self.other._getDiffusionGeomCoeff(mesh))

    def __neg__(self):
        r"""
         Negate a `_BinaryTerm`.

           >>> -(__ScalarCoeffTerm(coeff=1.) - __ScalarCoeffTerm(coeff=2.))
           (__ScalarCoeffTerm(coeff=-1.0) + __ScalarCoeffTerm(coeff=2.0))

        """

        return (-self.term) + (-self.other)

from fipy.terms.scalarCoeffTerm import _ScalarCoeffTerm
class __ScalarCoeffTerm(_ScalarCoeffTerm):
    """
    Dummy subclass for tests
    """
    pass 

def _test(): 
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
