#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "term.py"
 #                                    created: 11/12/03 {10:54:37 AM} 
 #                                last update: 10/26/06 {2:05:21 PM} 
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #  
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this software is not subject to copyright
 # protection and is in the public domain.  FiPy is an experimental
 # system.  NIST assumes no responsibility whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 # 
 # This software can be redistributed and/or modified freely
 # provided that any derivative works bear some notice that they are
 # derived from it, and any modified versions bear some notice that
 # they have been modified.
 # ========================================================================
 #  
 #  Description: 
 # 
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2003-11-12 JEG 1.0 original
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

import Numeric
from fipy.tools import numerix

from fipy.variables.variable import Variable
from fipy.tools.dimensions.physicalField import PhysicalField

class Term:
    """
    .. attention:: This class is abstract. Always create one of its subclasses.
    """
    def __init__(self, coeff=1.):
        """
        Create a `Term`.

        :Parameters:
          - `coeff`: The coefficient for the term. A `CellVariable` or number.
            `FaceVariable` objects are also acceptable for diffusion or convection terms.

        """  
        self.coeff = coeff
	self.geomCoeff = None
        self._cacheMatrix = False
        self.matrix = None
        self._cacheRHSvector = False
        self.RHSvector = None
        
    def _buildMatrix(self, var, boundaryConditions, dt):
	pass

    def _calcResidual(self, var, matrix, RHSvector):

	Lx = matrix * Numeric.array(var[:])
      
	residual = Lx - RHSvector
      
## 	denom = max(abs(Lx))
## 	if denom == 0:
## 	    denom = max(abs(RHSvector))
## 	if denom == 0:
## 	    denom = 1.
##           
## 	residual /= denom
      	
	return numerix.max(abs(residual))

    def __buildMatrix(self, var, boundaryConditions, dt):

        if type(boundaryConditions) not in (type(()), type([])):
            boundaryConditions = (boundaryConditions,)

        return self._buildMatrix(var, boundaryConditions, dt)

    def _solveLinearSystem(self, var, solver, matrix, RHSvector):

        if self._cacheMatrix:
            self.matrix = matrix

        if self._cacheRHSvector:
            self.RHSvector = RHSvector

        from fipy.solvers.linearPCGSolver import LinearPCGSolver

        solver = self._getDefaultSolver(solver) or solver or LinearPCGSolver()
	array = var.getNumericValue()
	solver._solve(matrix, array, RHSvector)
	var[:] = array
    
    def solve(self, var, solver=None, boundaryConditions=(), dt=1.):
        r"""
        Builds and solves the `Term`'s linear system once. This method
        does not return the residual. It should be used when the
        residual is not required.
      	
        :Parameters:

           - `var`: The variable to be solved for. Provides the initial condition, the old value and holds the solution on completion.
           - `solver`: The iterative solver to be used to solve the linear system of equations. Defaults to `LinearPCGSolver`.
           - `boundaryConditions`: A tuple of boundaryConditions.
           - `dt`: The time step size.

	"""
        self._verifyCoeffType(var)
        
        matrix, RHSvector = self.__buildMatrix(var, boundaryConditions, dt)
        
        self._solveLinearSystem(var, solver, matrix, RHSvector)

    def sweep(self, var, solver = None, boundaryConditions=(), dt=1., underRelaxation=None, residualFn = None):
        r"""
        Builds and solves the `Term`'s linear system once. This method
        also recalculates and returns the residual as well as applying
        under-relaxation.

        :Parameters:

           - `var`: The variable to be solved for. Provides the initial condition, the old value and holds the solution on completion.
           - `solver`: The iterative solver to be used to solve the linear system of equations. Defaults to `LinearPCGSolver`.
           - `boundaryConditions`: A tuple of boundaryConditions.
           - `dt`: The time step size.
           - `underRelaxation`: Usually a value between `0` and `1` or `None` in the case of no under-relaxation

	"""
        self._verifyCoeffType(var)
        
        matrix, RHSvector = self.__buildMatrix(var, boundaryConditions, dt)
 	
        if underRelaxation is not None:
            matrix, RHSvector = self._applyUnderRelaxation(matrix, var, RHSvector, underRelaxation)

        residualFn = residualFn or self._calcResidual
        residual = residualFn(var, matrix, RHSvector)
##         residual = self._calcResidual(var, matrix, RHSvector)

##         print "x", var

        self._solveLinearSystem(var, solver, matrix, RHSvector)

##         print "L", matrix
##         print "b", RHSvector

        return residual

    def _verifyCoeffType(self, var):
        pass

    def cacheMatrix(self):
        r"""
        Informs `solve()` and `sweep()` to cache their matrix so
        that `getMatrix()` can return the matrix.
        """
        self._cacheMatrix = True

    def getMatrix(self):
        r"""
        Return the matrix caculated in `solve()` or `sweep()`. The
        cacheMatrix() method should be called before `solve()` or
        `sweep()` to cache the matrix.
        """
        if not self._cacheMatrix: 
            import warnings
            warnings.warn("""cacheMatrix() should be called followed by sweep() or solve()
                          to calculate the matrix. None will be returned on this occasion.""",
                          UserWarning, stacklevel=2)
            self.cacheMatrix()

        return self.matrix

    def cacheRHSvector(self):
        r"""        
        Informs `solve()` and `sweep()` to cache their right hand side
        vector so that `getRHSvector()` can return it.
        """
        self._cacheRHSvector = True

    def getRHSvector(self):
        r"""
        Return the RHS vector caculated in `solve()` or `sweep()`. The
        cacheRHSvector() method should be called before `solve()` or
        `sweep()` to cache the vector.
        """
        if not self._cacheRHSvector: 
            import warnings
            warnings.warn("""getRHSvector should be called after cacheRHSvector() and either sweep() or solve()
                          to calculate the RHS vector. None will be returned on this occasion.""",
                          UserWarning, stacklevel=2)
            self.cacheRHSvector()

        return self.RHSvector
    


    def _applyUnderRelaxation(self, matrix, var, RHSVector, underRelaxation):
        matrix.putDiagonal(matrix.takeDiagonal() / underRelaxation)
        RHSVector += (1 - underRelaxation) * matrix.takeDiagonal() * numerix.array(var)
        return matrix, RHSVector

    def _getDefaultSolver(self, solver):
        return None
        
    def _otherIsZero(self, other):
        if (type(other) is type(0) or type(other) is type(0.)) and other == 0:
            return True
        else:
            return False

    def __add__(self, other):
        r"""
        Add a `Term` to another `Term`, number or variable.

           >>> Term(coeff = 1.) + 10.
           (Term(coeff = 1.0) + _ExplicitSourceTerm(coeff = 10.0))
           >>> Term(coeff = 1.) + Term(coeff = 2.)
           (Term(coeff = 1.0) + Term(coeff = 2.0))

        """
        
        if self._otherIsZero(other):
            return self
        else:
            from fipy.terms.binaryTerm import _AdditionTerm
            return _AdditionTerm(term1 = self, term2 = other)
	    
    __radd__ = __add__
    
    def __neg__(self):
        r"""
         Negate a `Term`.

           >>> -Term(coeff = 1.)
           Term(coeff = -1.0)

        """
        return self.__class__(coeff = -self.coeff)

    def __pos__(self):
        r"""
        Posate a `Term`.

           >>> +Term(coeff = 1.)
           Term(coeff = 1.0)

        """
        return self
        
    def __sub__(self, other):
        r"""
        Subtract a `Term` from a `Term`, number or variable.

           >>> Term(coeff = 1.) - 10.
           (Term(coeff = 1.0) - _ExplicitSourceTerm(coeff = 10.0))
           >>> Term(coeff = 1.) - Term(coeff = 2.)
           (Term(coeff = 1.0) - Term(coeff = 2.0))
           
        """        
        if self._otherIsZero(other):
            return self
        else:
            from fipy.terms.binaryTerm import _SubtractionTerm
            return _SubtractionTerm(term1 = self, term2 = other)

    def __rsub__(self, other):
        r"""
        Subtract a `Term`, number or variable from a `Term`.

           >>> 10. - Term(coeff = 1.)
           (_ExplicitSourceTerm(coeff = 10.0) - Term(coeff = 1.0))

        """        
        if self._otherIsZero(other):
            return self
        else:
            from fipy.terms.binaryTerm import _SubtractionTerm
            return _SubtractionTerm(term1 = other, term2 = self)
        
    def __eq__(self, other):
        r"""
        This method allows `Terms` to be equated in a natural way. Note that the
        following does not return `False.`

           >>> Term(coeff = 1.) == Term(coeff = 2.)
           (Term(coeff = 1.0) == Term(coeff = 2.0))

        it is equivalent to,

           >>> Term(coeff = 1.) - Term(coeff = 2.)
           (Term(coeff = 1.0) - Term(coeff = 2.0))

        A `Term` should equate with a float. 
        
        .. attention:: 
            
           This does not work due to sign difficulties.
           
        ..

           >>> Term(coeff = 1.) == 1.  
           False
           
        Likewise for integers.

           >>> Term(coeff = 1.) == 1
           False
           
        """

        if self._otherIsZero(other):
            return self
        else:
            if not isinstance(other, Term):
                return False
            else:
                from fipy.terms.binaryTerm import _EquationTerm
                return _EquationTerm(term1 = self, term2 = other)

                # because of the semantics of comparisons in Python,
                # the following test doesn't work
                ##         if isinstance(self, _EquationTerm) or isinstance(other, _EquationTerm):
                ##             raise SyntaxError, "Can't equate an equation with a term: %s == %s" % (str(self), str(other))

    def __repr__(self):
        """
        The representation of a `Term` object is given by,
        
           >>> print Term(123.456)
           Term(coeff = 123.456)

        """
        return "%s(coeff = %s)" % (self.__class__.__name__, str(self.coeff))

    def _calcGeomCoeff(self, mesh):
	return None
	
    def _getGeomCoeff(self, mesh):
	if self.geomCoeff is None:
	    self.geomCoeff = self._calcGeomCoeff(mesh)
            if self.geomCoeff is not None:
                self.geomCoeff.dontCacheMe()

	return self.geomCoeff
	
    def _getWeight(self, mesh):
	pass
	    
def _test(): 
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
