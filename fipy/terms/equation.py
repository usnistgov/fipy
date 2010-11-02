#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - a finite volume PDE solver in Python
 # 
 #  FILE: "equation.py"
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
from fipy.terms import TransientTerm, DiffusionTerm, \
  ExplicitDiffusionTerm, ImplicitSourceTerm
from fipy.terms.convectionTerm import ConvectionTerm
from fipy.terms.explicitSourceTerm import _ExplicitSourceTerm

class _Equation(Term):
    def __init__(self):
        self.orderedKeys = ["TransientTerm", "ExplicitDiffusionTerm", 
          "DiffusionTerm", "ConvectionTerm", "ImplicitSourceTerm", 
          "_ExplicitSourceTerm"]

        self.vars = []
        self.terms = []
        self.nonAdditiveTerms = []

	Term.__init__(self)
	
    def copy(self):
        eq = _Equation()
        eq.vars = self.vars.copy()
        for terms in self.terms:
            dup = {}
            for key, term in terms.iteritems():
                dup[key] = term.copy()
            eq.terms.append(dup)
        for terms in self.nonAdditiveTerms:
            dup = []
            for term in terms:
                dup.append(term.copy())
            eq.nonAdditiveTerms.append(dup)
        return eq
        
    def _getTermList(self, terms, nonAdditiveTerms):
        orderedPlusOtherKeys = self.orderedKeys \
          + [k for k in terms.keys() if k not in self.orderedKeys]
        return [terms[key] for key in orderedPlusOtherKeys if terms.has_key(key)] + nonAdditiveTerms
        
    def _getDiffusiveGeomCoeff(self, mesh):
        if self.terms.has_key("DiffusionTerm"):
            return self.terms["DiffusionTerm"]._getGeomCoeff(mesh)[0]
        else:
            return None
        
    def _checkAndBuildMatrix(self, var, SparseMatrix,  boundaryConditions, dt, equation=None):
        from fipy.tools import numerix

        self.RHSvector = 0
        self.matrix = 0

        for var, terms, nonAdditiveTerms in zip(self.vars, self.terms, self.nonAdditiveTerms):
            for term in self._getTermList(terms, nonAdditiveTerms):
                if term is not None:
                    var, termMatrix, termRHSvector = term._checkAndBuildMatrix(var, SparseMatrix,
                                                                               boundaryConditions, 
                                                                               dt, self)
                    
                    if (os.environ.has_key('FIPY_DISPLAY_MATRIX') 
                        and os.environ['FIPY_DISPLAY_MATRIX'].lower() == "terms"):
                        self._viewer.title = "%s %s" % (var.name, term.__class__.__name__)
                        self._viewer.plot(matrix=termMatrix, RHSvector=termRHSvector)
                        raw_input()

                    self.matrix += termMatrix
                    self.RHSvector += termRHSvector
                
        matrix = self.matrix
        RHSvector = self.RHSvector
        if not self._cacheMatrix:
            self.matrix = None
        if not self._cacheRHSvector:
            self.RHSvector = None
            
	return (var, matrix, RHSvector)
        
    def _getDefaultSolver(self, solver, *args, **kwargs):
        for var, terms, nonAdditiveTerms in zip(self.vars, self.terms, self.nonAdditiveTerms):
            for term in self._getTermList(terms, nonAdditiveTerms):
                defaultsolver = term._getDefaultSolver(solver, *args, **kwargs)
                if defaultsolver is not None:
                    return defaultsolver
                
        return solver

    def __repr__(self):
        reprs = []
        for var, terms, nonAdditiveTerms in zip(self.vars, self.terms, self.nonAdditiveTerms):
            for term in self._getTermList(terms, nonAdditiveTerms):
                reprs.append(repr(term))
                
        return " + ".join(reprs) + " == 0"

    def __add__(self, other):
        r"""
        Add a `Term` to another `Term`, number or variable.

           >>> __Term(coeff=1.) + 10. + __Term(2.)
           10.0 + __Term(coeff=3.0) == 0
           >>> __Term(coeff=1.) + __Term(coeff=2.) + __Term(coeff=3.)
           __Term(coeff=6.0)

        """
        if self._otherIsZero(other):
            return self
        else:
            dup = self.copy()
            dup += other
                
            return dup
        
    def _appendTerm(self, other, key):
        if self.terms[key] is None:
            self.terms[key] = other
        else:
            self.terms[key] = self.terms[key] + other
        
    def __iadd__(self, other):
        if other is not None:
            if not isinstance(other, Term):
                self._appendTerm(_ExplicitSourceTerm(coeff=other), "_ExplicitSourceTerm")
            elif isinstance(other, _Equation):
                for terms in other.terms:
                    for key in terms.keys():
                        self += terms[key]
            else:
                appended = False

                try:
                    # find-in-list generator: http://dev.ionous.net/2009/01/python-find-item-in-list.html
                    var_index = (i for i, var in enumerate(self.vars) if var is other.var).next()
                except StopIteration:
                    if other.var is None and len(self.vars) > 0:
                        raise Exception("Terms with explicit Variables cannot mix with Terms with implicit Variables")

                    self.vars.append(other.var)
                    self.terms.append({})
                    self.nonAdditiveTerms.append([])
                    var_index = len(self.vars) - 1
                
                for key in self.terms[var_index].keys():

                    if eval("isinstance(other, %s)" % key):
                        self._appendTerm(other, key)
                        appended = True
                        break
                        
                if not appended:
                    if other._isAdditive():
                        self.terms[var_index][other.__class__.__name__] = other
                    else:

                        self.nonAdditiveTerms[var_index] += [other]
        
        return self
        
    def __neg__(self):
        r"""
         Negate a `Term`.

           >>> -(__Term(coeff=1.) - __Term(coeff=2.))
           __Term(coeff=1.0)

        """
        dup = self.copy()
        
        for key, term in self.terms.iteritems():
            if term is not None:
                dup.terms[key] = -term
                
        return dup

    def _test(self):
        """
        These tests are not useful as documentation, but are here to ensure
        everything works as expected.
        
        >>> from fipy import Grid1D, CellVariable, DiffusionTerm, TransientTerm
        >>> mesh = Grid1D(nx=3)
        >>> A = CellVariable(mesh=mesh, name="A")
        >>> B = CellVariable(mesh=mesh, name="B")
        >>> eq = TransientTerm(coeff=1., var=A) == DiffusionTerm(coeff=1., var=B)
        >>> print eq
        TransientTerm(coeff=1.0, var=A) + DiffusionTerm(coeff=[-1.0], var=B) == 0
        >>> print eq.vars
        >>> print eq.terms
        >>> solver = eq._prepareLinearSystem(var=None, solver=None, boundaryConditions=(), dt=1.)
        >>> print solver.matrix
        >>> print solver.RHSvector
        """

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

