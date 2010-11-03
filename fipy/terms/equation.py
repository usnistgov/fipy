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
        eq.vars = []
        eq.vars.extend(self.vars)
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
        
    def _findVarIndex(self, var):
        # find-in-list generator: http://dev.ionous.net/2009/01/python-find-item-in-list.html
        return (i for i, v in enumerate(self.vars) if v is var).next()
        
    def _findOrAppendVarIndex(self, var):
        try:
            var_index = self._findVarIndex(var=var)
        except StopIteration:
            self.vars.append(var)
            self.terms.append({})
            self.nonAdditiveTerms.append([])
            var_index = len(self.vars) - 1
            
        return var_index

    def _getOrderdPlusOtherKeys(self, terms):
        return self.orderedKeys + [k for k in terms.keys() if k not in self.orderedKeys]

    def _getTermList(self, terms, nonAdditiveTerms):
        orderedPlusOtherKeys = self._getOrderdPlusOtherKeys(terms)
        return [terms[key] for key in orderedPlusOtherKeys if terms.has_key(key)] + nonAdditiveTerms
        
    def _checkAndBuildMatrix(self, var, SparseMatrix,  boundaryConditions, dt, equation=None):
        from fipy.tools import numerix

        self.RHSvector = 0
        self.matrix = 0
        
        try:
            var_index = self._findVarIndex(var=var)
        except StopIteration:
            if len(self.vars) == 0:
                # equation is not defined on any Variables
                return (var, None, None)
            elif len(self.vars) == 1:
                if var is None:
                    var_index = 0
                else:
                    var_index = None
            else:
                if var is None:
                    raise Exception("Can't build matrix without specifying a Variable")
                else:
                    # equation is not defined on this Variable
                    return (var, None, None)

        if var_index is None:
            var_index = 0
        else:
            var = self.vars[var_index]
        
        if var is None:
            raise Exception("Can't build matrix without specifying a Variable")

        terms = self.terms[var_index]
        nonAdditiveTerms = self.nonAdditiveTerms[var_index]
        for term in self._getTermList(terms, nonAdditiveTerms):
            if self.matrix == 0:
                term._setDiagonalSign(sign=1)
            else:
                from fipy.tools.numerix import sign, add
                term._setDiagonalSign(sign=sign(add.reduce(self.matrix.takeDiagonal())))
                
            if terms.has_key("DiffusionTerm"):
                term._setDiffusiveGeomCoeff(diffCoeff=terms["DiffusionTerm"]._getGeomCoeff(var.getMesh())[0])
            else:
                term._setDiffusiveGeomCoeff(diffCoeff=None)

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
        var_index = self._findOrAppendVarIndex(var=other.var)

        terms = self.terms[var_index]
        if terms.has_key(key):
            terms[key] = terms[key] + other
        else:
            terms[key] = other
        
    def __iadd__(self, other):
        if not isinstance(other, Term):
            self._appendTerm(_ExplicitSourceTerm(coeff=other), "_ExplicitSourceTerm")
        elif isinstance(other, _Equation):
            for terms in other.terms:
                for key in terms.keys():
                    self += terms[key]
        else:
            appended = False

            if other.var is None and len(self.vars) > 0:
                if ((len(self.vars) == 1 and self.vars[0] is not None)
                    or (len(self.vars) > 1 and not isinstance(other, _ExplicitSourceTerm))):
                        raise Exception("Terms with explicit Variables cannot mix with Terms with implicit Variables")

            var_index = self._findOrAppendVarIndex(var=other.var)
            
            for key in self._getOrderdPlusOtherKeys(self.terms[var_index]):
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
        
        for terms in dup.terms:
            for key, term in terms.iteritems():
                terms[key] = -term
            
        dup.nonAdditiveTerms = [[-term for term in terms] for terms in dup.nonAdditiveTerms]
                
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
        [A, B]
        >>> print eq.terms
        [{'TransientTerm': TransientTerm(coeff=1.0, var=A)}, {'DiffusionTerm': DiffusionTerm(coeff=[-1.0], var=B)}]
        >>> solver = eq._prepareLinearSystem(var=None, solver=None, boundaryConditions=(), dt=1.)
        Traceback (most recent call last):
            ...
        Exception: Can't build matrix without specifying a Variable
        >>> solver = eq._prepareLinearSystem(var=A, solver=None, boundaryConditions=(), dt=1.)
        >>> print solver.matrix
         1.000000      ---        ---    
            ---     1.000000      ---    
            ---        ---     1.000000  
        >>> print solver.RHSvector
        [ 0.  0.  0.]
        >>> solver = eq._prepareLinearSystem(var=B, solver=None, boundaryConditions=(), dt=1.)
        >>> print solver.matrix
         1.000000  -1.000000      ---    
        -1.000000   2.000000  -1.000000  
            ---    -1.000000   1.000000  
        >>> print solver.RHSvector
        [ 0.  0.  0.]
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

