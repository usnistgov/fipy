#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - a finite volume PDE solver in Python
 # 
 #  FILE: "binaryTerm.py"
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
from fipy.terms.baseBinaryTerm import _BaseBinaryTerm
from fipy.terms.explicitSourceTerm import _ExplicitSourceTerm

class _BinaryTerm(_BaseBinaryTerm):

    def _verifyVar(self, var):

        if var is None and len(self._getVars()) > 1:
            raise Exception, 'The solution variable needs to be specified'

        return Term._verifyVar(self, var)
    
    def _buildMatrix(self, var, SparseMatrix,  boundaryConditions=(), dt=1.0, transientGeomCoeff=None, diffusionGeomCoeff=None):

        matrix = 0
        RHSvector = 0

        for term in (self.term, self.other):

            termVar, termMatrix, termRHSvector = term._buildMatrix(var,
                                                                   SparseMatrix,
                                                                   boundaryConditions=boundaryConditions,
                                                                   dt=dt,
                                                                   transientGeomCoeff=transientGeomCoeff,
                                                                   diffusionGeomCoeff=diffusionGeomCoeff)

            matrix += termMatrix
            RHSvector += termRHSvector

        if (os.environ.has_key('FIPY_DISPLAY_MATRIX')
            and os.environ['FIPY_DISPLAY_MATRIX'].lower() == "terms"): 
            self._viewer.title = "%s %s" % (var.name, repr(self))
            self._viewer.plot(matrix=matrix, RHSvector=RHSvector) 
            raw_input()

	return (var, matrix, RHSvector)
        
    def __repr__(self):
        return '(' + repr(self.term) + ' + ' + repr(self.other) + ')'

    def __neg__(self):
        r"""
         Negate a `_BinaryTerm`.

           >>> -(__Term(coeff=1.) - __Term(coeff=2.))
           (__Term(coeff=-1.0) + __Term(coeff=2.0))

        """

        return (-self.term) + (-self.other)

    def __mul__(self, other):
        return other * self.term + other * self.other

    __rmul__ = __mul__

class __Term(Term):
    """
    Dummy subclass for tests
    """
    pass 

def _test(): 
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
