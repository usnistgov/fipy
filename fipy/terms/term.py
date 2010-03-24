#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "term.py"
 #
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
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

import os

from fipy.tools import numerix

from fipy.variables.variable import Variable
from fipy.tools.dimensions.physicalField import PhysicalField
from fipy.solvers import DefaultSolver

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
        if self.__class__ is Term:
            raise NotImplementedError, "can't instantiate abstract base class"
            
        self.coeff = coeff
        self.geomCoeff = None
        self._cacheMatrix = False
        self.matrix = None
        self._cacheRHSvector = False
        self.RHSvector = None
        self._diagonalSign = Variable(value=1)
        
    def copy(self):
        return self.__class__(self.coeff)
        
    def _buildMatrix(self, var, SparseMatrix, boundaryConditions, dt, equation=None):
        raise NotImplementedError

    def __buildMatrix(self, var, solver, boundaryConditions, dt):
        if numerix.sctype2char(var.getsctype()) not in numerix.typecodes['Float']:
            import warnings
            warnings.warn("""sweep() or solve() are likely to produce erroneous results when `var` does not contain floats.""",
                          UserWarning, stacklevel=4)
        
        self._verifyCoeffType(var)
        
        if numerix.getShape(dt) != ():
            raise TypeError, "`dt` must be a single number, not a " + type(dt).__name__
        dt = float(dt)
    
        if type(boundaryConditions) not in (type(()), type([])):
            boundaryConditions = (boundaryConditions,)

        for bc in boundaryConditions:
            bc._resetBoundaryConditionApplied()

        if os.environ.has_key('FIPY_DISPLAY_MATRIX'):
            if not hasattr(self, "_viewer"):
                from fipy.viewers.matplotlibViewer.matplotlibSparseMatrixViewer import MatplotlibSparseMatrixViewer
                Term._viewer = MatplotlibSparseMatrixViewer()

        matrix, RHSvector = self._buildMatrix(var, solver._getMatrixClass(), boundaryConditions, dt)
        
        solver._storeMatrix(var=var, matrix=matrix, RHSvector=RHSvector)
        
        if os.environ.has_key('FIPY_DISPLAY_MATRIX'):
            self._viewer.title = "%s %s" % (var.name, self.__class__.__name__)
            self._viewer.plot(matrix=matrix, RHSvector=RHSvector)
            from fipy import raw_input
            raw_input()

    def _prepareLinearSystem(self, var, solver, boundaryConditions, dt):
        solver = self.getDefaultSolver(solver)

        self.__buildMatrix(var, solver, boundaryConditions, dt)
        return solver
    
    def solve(self, var, solver=None, boundaryConditions=(), dt=1.):
        r"""
        Builds and solves the `Term`'s linear system once. This method
        does not return the residual. It should be used when the
        residual is not required.
        
        :Parameters:

           - `var`: The variable to be solved for. Provides the initial condition, the old value and holds the solution on completion.
           - `solver`: The iterative solver to be used to solve the linear system of equations. Defaults to `LinearPCGSolver` for Pysparse and `LinearLUSolver` for Trilinos.
           - `boundaryConditions`: A tuple of boundaryConditions.
           - `dt`: The time step size.

        """
        
        solver = self._prepareLinearSystem(var, solver, boundaryConditions, dt)
        
        solver._solve()

    def sweep(self, var, solver = None, boundaryConditions=(), dt=1., underRelaxation=None, residualFn=None):
        r"""
        Builds and solves the `Term`'s linear system once. This method
        also recalculates and returns the residual as well as applying
        under-relaxation.

        :Parameters:

           - `var`: The variable to be solved for. Provides the initial condition, the old value and holds the solution on completion.
           - `solver`: The iterative solver to be used to solve the linear system of equations. Defaults to `LinearPCGSolver` for Pysparse and `LinearLUSolver` for Trilinos.
           - `boundaryConditions`: A tuple of boundaryConditions.
           - `dt`: The time step size.
           - `underRelaxation`: Usually a value between `0` and `1` or `None` in the case of no under-relaxation
           - `residualFn`: A function that takes var, matrix, and RHSvector arguments, used to customize the residual calculation.

        """
        solver = self._prepareLinearSystem(var=var, solver=solver, boundaryConditions=boundaryConditions, dt=dt)
        solver._applyUnderRelaxation(underRelaxation=underRelaxation)
        residual = solver._calcResidual(residualFn=residualFn)

        solver._solve()

        return residual

    def justResidualVector(self, var, solver=None, boundaryConditions=(), dt=1., underRelaxation=None, residualFn=None):
        r"""
        Builds the `Term`'s linear system once. This method
        also recalculates and returns the residual as well as applying
        under-relaxation.

        :Parameters:

           - `var`: The variable to be solved for. Provides the initial condition, the old value and holds the solution on completion.
           - `solver`: The iterative solver to be used to solve the linear system of equations. Defaults to `LinearPCGSolver` for Pysparse and `LinearLUSolver` for Trilinos.
           - `boundaryConditions`: A tuple of boundaryConditions.
           - `dt`: The time step size.
           - `underRelaxation`: Usually a value between `0` and `1` or `None` in the case of no under-relaxation
           - `residualFn`: A function that takes var, matrix, and RHSvector arguments used to customize the residual calculation.

        """
        solver = self._prepareLinearSystem(var, solver, boundaryConditions, dt)
        solver._applyUnderRelaxation(underRelaxation)

        return solver._calcResidualVector(residualFn=residualFn)

    def residualVectorAndNorm(self, var, solver=None, boundaryConditions=(), dt=1., underRelaxation=None, residualFn=None):
        r"""
        Builds the `Term`'s linear system once. This method
        also recalculates and returns the residual as well as applying
        under-relaxation.

        :Parameters:

           - `var`: The variable to be solved for. Provides the initial condition, the old value and holds the solution on completion.
           - `solver`: The iterative solver to be used to solve the linear system of equations. Defaults to `LinearPCGSolver` for Pysparse and `LinearLUSolver` for Trilinos.
           - `boundaryConditions`: A tuple of boundaryConditions.
           - `dt`: The time step size.
           - `underRelaxation`: Usually a value between `0` and `1` or `None` in the case of no under-relaxation
           - `residualFn`: A function that takes var, matrix, and RHSvector arguments used to customize the residual calculation.

        """
        solver = self._prepareLinearSystem(var, solver, boundaryConditions, dt)
        solver._applyUnderRelaxation(underRelaxation)
        vector = solver._calcResidualVector(residualFn=residualFn)
        
        L2norm = numerix.L2norm(vector)

        return vector, L2norm

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
    
    def _getDefaultSolver(self, solver, *args, **kwargs):
        return None
        
    def getDefaultSolver(self, solver=None, *args, **kwargs):
        return self._getDefaultSolver(solver, *args, **kwargs) or solver or DefaultSolver(*args, **kwargs)
                         
    def _otherIsZero(self, other):
        if (type(other) is type(0) or type(other) is type(0.)) and other == 0:
            return True
        else:
            return False
            
    def __add__(self, other):
        r"""
        Add a `Term` to another `Term`, number or variable.

           >>> __Term(coeff=1.) + 10.
           10.0 + __Term(coeff=1.0) == 0
           >>> __Term(coeff=1.) + __Term(coeff=2.)
           __Term(coeff=3.0)

        """
        from fipy.terms.equation import _Equation

##        print 'self',self
##        print 'other',other
##        print 'isinstance(other, _Equation)',isinstance(other, _Equation)
        
        if self._otherIsZero(other):
            return self
        elif isinstance(other, _Equation):
            return other + self
        elif self.__class__ == other.__class__:
            return self.__class__(coeff=self.coeff + other.coeff)
        else:
            return self._add(other)
                
    def _add(self, other):
        from fipy.terms.equation import _Equation
        eq = _Equation()
        eq += self
        eq += other
        return eq
            
    def __radd__(self, other):
        r"""
        Add a number or variable to a `Term`.

           >>> 10. + __Term(coeff=1.)
           10.0 + __Term(coeff=1.0) == 0
        """
        return self + other
    
    def __neg__(self):
        r"""
         Negate a `Term`.

           >>> -__Term(coeff=1.)
           __Term(coeff=-1.0)

        """
        try:
            coeff = -self.coeff
        except:
            coeff = -numerix.array(self.coeff)

        return self.__class__(coeff=coeff)

    def __pos__(self):
        r"""
        Posate a `Term`.

           >>> +__Term(coeff=1.)
           __Term(coeff=1.0)

        """
        return self
        
    def __sub__(self, other):
        r"""
        Subtract a `Term` from a `Term`, number or variable.

           >>> __Term(coeff=1.) - 10.
           -10.0 + __Term(coeff=1.0) == 0
           >>> __Term(coeff=1.) - __Term(coeff=2.)
           __Term(coeff=-1.0)
           
        """        
        if self._otherIsZero(other):
            return self
        else:
            return self + (-other)

    def __rsub__(self, other):
        r"""
        Subtract a `Term`, number or variable from a `Term`.

           >>> 10. - __Term(coeff=1.)
           10.0 + __Term(coeff=-1.0) == 0

        """        
        if self._otherIsZero(other):
            return -self
        else:
            return other + (-self)
        
    def __eq__(self, other):
        r"""
        This method allows `Terms` to be equated in a natural way. Note that the
        following does not return `False.`

           >>> __Term(coeff=1.) == __Term(coeff=2.)
           __Term(coeff=-1.0)

        it is equivalent to,

           >>> __Term(coeff=1.) - __Term(coeff=2.)
           __Term(coeff=-1.0)

        A `Term` can also equate with a number. 

           >>> __Term(coeff=1.) == 1.  
           -1.0 + __Term(coeff=1.0) == 0
           
        Likewise for integers.

           >>> __Term(coeff=1.) == 1
           -1 + __Term(coeff=1.0) == 0
           
        Equating to zero is allowed, of course
        
            >>> __Term(coeff=1.) == 0
            __Term(coeff=1.0)
            >>> 0 == __Term(coeff=1.)
            __Term(coeff=1.0)
           
        """

        if self._otherIsZero(other):
            return self
        else:
            return self - other

    def __mul__(self, other):
        r"""
        Mutiply a term

            >>> 2. * __Term(coeff=0.5)
            2.0 * __Term(coeff=0.5)
            
        """         
        from fipy.terms.mulTerm import _MulTerm
        return _MulTerm(term=self, coeff=other)

    __rmul__ = __mul__
               
    def __div__(self, other):
        r"""
        Divide a term

            >>> __Term(2.) / 2.
            0.5 * __Term(coeff=2.0)

        """
        return (1 / other) * self
    
    def __repr__(self):
        """
        The representation of a `Term` object is given by,
        
           >>> print __Term(123.456)
           __Term(coeff=123.456)

        """
        return "%s(coeff=%s)" % (self.__class__.__name__, repr(self.coeff))

    def _calcGeomCoeff(self, mesh):
        return None
        
    def _getGeomCoeff(self, mesh):
        if self.geomCoeff is None:
            self.geomCoeff = self._calcGeomCoeff(mesh)
            if self.geomCoeff is not None:
                self.geomCoeff.dontCacheMe()

        return self.geomCoeff
        
    def _getWeight(self, mesh):
        raise NotImplementedError

    def _isAdditive(self):
        return True
            
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
