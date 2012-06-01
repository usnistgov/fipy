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

__docformat__ = 'restructuredtext'

__all__ = []

import os

from fipy.terms.term import Term
from fipy.terms.explicitSourceTerm import _ExplicitSourceTerm
from fipy.terms import ExplicitVariableError

class _AbstractBinaryTerm(Term):
    def __init__(self, term, other):

        if not isinstance(other, Term):
            other = _ExplicitSourceTerm(coeff=other, var=term.var)

        self.term = term
        self.other = other

        if term.var is None:
            if other.var is None:
                pass
            else:
                raise ExplicitVariableError
        else:
            if other.var is None:
                if isinstance(other, _ExplicitSourceTerm):
                    other.var = term.var
                else:
                    raise ExplicitVariableError

        Term.__init__(self, var=self._vars[0])
    
    def _addNone(self, arg0, arg1):
        if arg0 is None and arg1 is None:
            return None
        elif arg0 is None:
            return arg1
        elif arg1 is None:
            return arg0
        else:
            return arg0 + arg1

    def __neg__(self):
        r"""
         Negate a `_BinaryTerm`.

           >>> -(__NonDiffusionTerm(coeff=1.) - __NonDiffusionTerm(coeff=2.))
           (__NonDiffusionTerm(coeff=-1.0) + __NonDiffusionTerm(coeff=2.0))

        """

        return (-self.term) + (-self.other)

    def _calcVars(self):
        return list(set(self.term._vars + self.other._vars))

    @property
    def _vars(self):
        if not hasattr(self, '_internalVars'):
            self._internalVars = self._calcVars()
        return self._internalVars

    @property
    def _transientVars(self):
        return self.term._transientVars + self.other._transientVars

    @property
    def _diffusionVars(self):
        return self.term._diffusionVars + self.other._diffusionVars

    def _checkVar(self, var):
        self.term._checkVar(var)
        self.other._checkVar(var)
    

from fipy.terms.nonDiffusionTerm import _NonDiffusionTerm
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
