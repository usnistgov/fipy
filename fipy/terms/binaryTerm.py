#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - a finite volume PDE solver in Python
 # 
 #  FILE: "binaryTerm.py"
 #                                    created: 11/9/04 {11:51:08 AM} 
 #                                last update: 3/23/07 {7:49:25 AM} 
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
 #  Description: 
 # 
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2004-11-09 JEG 1.0 original
 # ###################################################################
 ##

from fipy.terms.term import Term
from fipy.terms.explicitSourceTerm import _ExplicitSourceTerm

class _AdditionTerm(Term):
    def __init__(self, term1, term2):

	if isinstance(term1, Term):
	    if not isinstance(term2, Term):
		term2 = _ExplicitSourceTerm(coeff = term2)
	elif isinstance(term2, Term):
	    term1 = _ExplicitSourceTerm(coeff = term1)
	else:
	    raise "No terms!"
	
	self.term1 = term1
	self.term2 = term2
	    
	Term.__init__(self)
	
    def _buildMatrix(self, var, boundaryConditions, dt):
	matrix, RHSvector = self.term1._buildMatrix(var, boundaryConditions, dt = dt)
	
	termMatrix, termRHSvector = self.term2._buildMatrix(var, boundaryConditions, dt = dt)

	matrix = matrix + termMatrix
	RHSvector = RHSvector + termRHSvector
	
	return (matrix, RHSvector)
        
    def _getDefaultSolver(self, solver):
        return self.term1._getDefaultSolver(solver) or self.term2._getDefaultSolver(solver)

    def __repr__(self):
        return "%s + %s" % (repr(self.term1), repr(self.term2))

    def __neg__(self):
        r"""
         Negate a `Term`.

           >>> -(Term(coeff=1.) - Term(coeff=2.)
           (Term(coeff = -1.0) + Term(coeff = 2.0))

        """
        return (-self.term1) + (-self.term2)

## class _AdditionTerm(_BinaryTerm):
##     def __repr__(self):
##         return "(%s + %s)" % (repr(self.term1), repr(self.term2))
##         
##     def _operator(self):
## 	return lambda a,b: a + b
## 	
## class _SubtractionTerm(_BinaryTerm):
##     def __repr__(self):
##         return "(%s - %s)" % (repr(self.term1), repr(self.term2))
##         
##     def _operator(self):
## 	return lambda a,b: a - b
## 
## class _EquationTerm(_SubtractionTerm):
##     def __repr__(self):
##         return "(%s == %s)" % (repr(self.term1), repr(self.term2))
##         
##     def __nonzero__(self):
##         if self.term1.__class__ != self.term2.__class__:
##             return False
##         elif self.term1.coeff != self.term2.coeff:
##             return False
##         else:
##             return True
##         
