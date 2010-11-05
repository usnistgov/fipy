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
from fipy.terms.explicitSourceTerm import _ExplicitSourceTerm

class _BinaryTerm(Term):
    def __init__(self, term, other):
        if not isinstance(other, Term):
            other = _ExplicitSourceTerm(other)
        self.terms = (term, other)
	Term.__init__(self)
	
    def _buildMatrix(self, var, SparseMatrix,  boundaryConditions=(), dt=1.0)

        matrix = 0
        RHSvector = 0

        for term in self.terms:
            tmpMatrix, tmpRHSvector = term._buildMatrix(var,
                                                        SparseMatrix,
                                                        boundaryConditions=boundaryConditions,
                                                        dt=dt)
##            from fipy.tools.debug import PRINT
##            PRINT('term')
##            PRINT('matrix',matrix)
##            PRINT('tmpMatrix',tmpMatrix)
            matrix += tmpMatrix
            RHSvector += tmpRHSvector
##            PRINT('matrix',matrix)
##            raw_input('stopped')
	return (matrix, RHSvector)
        
    def _getDefaultSolver(self, solver, *args, **kwargs):
         for term in self.terms:
             defaultsolver = term._getDefaultSolver(solver, *args, **kwargs)
             if defaultsolver is not None:
                 return defaultsolver
                
         return solver

    def __repr__(self):

        return '(' + repr(self.terms[0]) + ' + ' + repr(self.terms[1]) + ')'

    def __neg__(self):
        r"""
         Negate a `_BinaryTerm`.

           >>> -(__Term(coeff=1.) - __Term(coeff=2.))
           (__Term(coeff=-1.0) + __Term(coeff=2.0))

        """

        return (-self.terms[0]) + (-self.terms[1])

    def __mul__(self, other):
        return other * self.terms[0] + other * self.terms[1]

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

