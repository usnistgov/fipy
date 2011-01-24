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
from fipy.solvers import DefaultSolver


class Term:
    """
    .. attention:: This class is abstract. Always create one of its subclasses.
    """
    self.__isAbstractBaseClass
    def __init__(self, coeff=1., var=None):
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
        self.var = var

    def _getVars(self):
        raise NotImplementedError

    def _calcVars(self):
        raise NotImplementedError

    def _getTransientVars(self):
        raise NotImplementedError

    def _getDiffusionVars(self):
        raise NotImplementedError
                
    def copy(self):
        return self.__class__(self.coeff, var=self.var)
        
    def _buildMatrix(self, var, SparseMatrix, boundaryConditions=(), dt=1.0, transientGeomCoeff=None, diffusionGeomCoeff=None):
        raise NotImplementedError

    def _verifyVar(self, var):
        if var is None:
            if self.var is None:
                raise Exception, 'The solution variable needs to be specified'
            else:
                return self.var
        else:
            return var

    def _checkVar(self, var):
        if numerix.sctype2char(var.getsctype()) not in numerix.typecodes['Float']:
            import warnings
            warnings.warn("""sweep() or solve() are likely to produce erroneous results when `var` does not contain floats.""",
                          UserWarning, stacklevel=4)
        
    def __buildMatrix(self, var, solver, boundaryConditions, dt):

        var = self._verifyVar(var)
        self._checkVar(var)

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

        var, matrix, RHSvector = self._buildMatrix(var, solver._getMatrixClass(), boundaryConditions, dt,
                                                   transientGeomCoeff=self._getTransientGeomCoeff(var),
                                                   diffusionGeomCoeff=self._getDiffusionGeomCoeff(var))
        
        if self._cacheMatrix:
            self.matrix = matrix
            self.matrix.cache = True
        else:
            self.matrix = None

        if self._cacheRHSvector:
            self.RHSvector = RHSvector
        else:
            self.RHSvector = None
        
        solver._storeMatrix(var=var, matrix=matrix, RHSvector=RHSvector)
        
        if (os.environ.has_key('FIPY_DISPLAY_MATRIX')
            and os.environ['FIPY_DISPLAY_MATRIX'].lower() == "terms"): 
            self._viewer.title = r"%s %s" % (var.name, self.__class__.__name__)
            self._viewer.plot(matrix=matrix, RHSvector=RHSvector)
            from fipy import raw_input
            raw_input()

    def _prepareLinearSystem(self, var, solver, boundaryConditions, dt):
        solver = self.getDefaultSolver(solver)
            
        self.__buildMatrix(var, solver, boundaryConditions, dt)
        
        if os.environ.has_key('FIPY_DISPLAY_MATRIX'):
            if var is None:
                name = ""
            else:
                name = var.name
            self._viewer.title = r"%s %s" % (name, repr(self))
            from fipy.variables.coupledCellVariable import _CoupledCellVariable
            if isinstance(solver.RHSvector, _CoupledCellVariable):
                RHSvector = solver.RHSvector.getGlobalValue()
            else:
                RHSvector = solver.RHSvector
            self._viewer.plot(matrix=solver.matrix, RHSvector=RHSvector)
            from fipy import raw_input
            raw_input()

        return solver
    
    def solve(self, var=None, solver=None, boundaryConditions=(), dt=1.):
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

    def sweep(self, var=None, solver=None, boundaryConditions=(), dt=1., underRelaxation=None, residualFn=None):
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

    def justResidualVector(self, var=None, solver=None, boundaryConditions=(), dt=1., underRelaxation=None, residualFn=None):
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

    def residualVectorAndNorm(self, var=None, solver=None, boundaryConditions=(), dt=1., underRelaxation=None, residualFn=None):
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
                         
    def __add__(self, other):
        if isinstance(other, (int, float)) and other == 0:
            return self
        else:
            from fipy.terms.binaryTerm import _BinaryTerm
            return _BinaryTerm(self, other)

    __radd__ = __add__
    
    def __neg__(self):
        raise NotImplementedError

    def __pos__(self):
        return self
        
    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return other + (-self)
        
    def __eq__(self, other):
        return self - other

    def __mul__(self, other):
        raise NotImplementedError

    __rmul__ = __mul__

    def __div__(self, other):
        return (1 / other) * self

    def __and__(self, other):
        if isinstance(other, Term):
            from fipy.terms.coupledBinaryTerm import _CoupledBinaryTerm
            return _CoupledBinaryTerm(self, other)
        else:
            raise Exception

    def _getUncoupledTerms(self):
        raise NotImplementedError
    
    def __repr__(self):
        raise NotImplementedError

    def _calcGeomCoeff(self, mesh):
        raise NotImplementedError
        
    def _getGeomCoeff(self, mesh):
        if self.geomCoeff is None:
            self.geomCoeff = self._calcGeomCoeff(mesh)
            if self.geomCoeff is not None:
                self.geomCoeff.dontCacheMe()

        return self.geomCoeff
        
    def _getWeight(self, var, transientGeomCoeff=None, diffusionGeomCoeff=None):
        raise NotImplementedError

    def _getDiffusionGeomCoeff(self, var):
        return None

    def _getTransientGeomCoeff(self, var):
        return None

    def _treatMeshAsOrthogonal(self, mesh):
        raise NotImplementedError

    def _getNormals(self, mesh):
        raise NotImplementedError

    def _getOldAdjacentValues(self, oldArray, id1, id2, dt):
        raise NotImplementedError

    def _getDifferences(self, adjacentValues, cellValues, oldArray, cellToCellIDs, mesh):
        raise NotImplementedError

    def _test(self):
        """
        Test stuff.
    
        >>> from fipy import *
        >>> L = 1.
        >>> nx = 100
        >>> m = Grid1D(nx=nx, dx=L / nx)
        >>> v = CellVariable(mesh=m, value=1.)
        >>> eqn = DiffusionTerm() - v
        
        >>> v.constrain(0.,  m.getFacesLeft())
        >>> v.constrain(1.,  m.getFacesRight())
        
        >>> res = 1.
        >>> sweep = 0
        >>> while res > 1e-8 and sweep < 100:
        ...     res = eqn.sweep(v)
        ...     sweep += 1
        >>> x = m.getCellCenters()[0]
        >>> answer = (numerix.exp(x) - numerix.exp(-x)) / (numerix.exp(L) - numerix.exp(-L))
        >>> print numerix.allclose(v, answer, rtol=2e-5)
        True
        
        >>> v.setValue(0.)
        >>> eqn = DiffusionTerm(0.2) * 5. - 5. * ImplicitSourceTerm(0.2)
        >>> eqn.solve(v)
        >>> print numerix.allclose(v, answer, rtol=2e-5)
        True
        
        >>> v.setValue(0.)
        >>> eqn = 2. * (DiffusionTerm(1.) - ImplicitSourceTerm(.5)) - DiffusionTerm(1.)
        >>> eqn.solve(v)
        >>> print numerix.allclose(v, answer, rtol=2e-5)
        True

 	>>> from fipy import Grid1D, CellVariable, DiffusionTerm, TransientTerm 
 	>>> mesh = Grid1D(nx=3) 
 	>>> A = CellVariable(mesh=mesh, name="A") 
 	>>> B = CellVariable(mesh=mesh, name="B") 
 	>>> eq = TransientTerm(coeff=1., var=A) == DiffusionTerm(coeff=1., var=B) 
 	>>> print eq 
 	(TransientTerm(coeff=1.0, var=A) + DiffusionTerm(coeff=[-1.0], var=B))
 	>>> A in set(eq._getVars()) and B in set(eq._getVars()) ## _getVars() is unordered for _BinaryTerm's.
 	True
 	>>> print (eq.term, eq.other) 
 	(TransientTerm(coeff=1.0, var=A), DiffusionTerm(coeff=[-1.0], var=B))
 	>>> solver = eq._prepareLinearSystem(var=None, solver=None, boundaryConditions=(), dt=1.)
 	Traceback (most recent call last): 
 	    ... 
        Exception: The solution variable needs to be specified
 	>>> solver = eq._prepareLinearSystem(var=A, solver=None, boundaryConditions=(), dt=1.)
        >>> from fipy.tools import parallel
        >>> numpyMatrix = solver.matrix.getNumpyArray()
 	>>> print parallel.procID > 0 or numerix.allequal(numpyMatrix, [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        True
        >>> print parallel.procID > 0 or numerix.allequal(solver.RHSvector, [0, 0, 0])
 	True
 	>>> solver = eq._prepareLinearSystem(var=B, solver=None, boundaryConditions=(), dt=1.)
        >>> numpyMatrix = solver.matrix.getNumpyArray()
        >>> print parallel.procID > 0 or numerix.allequal(numpyMatrix, [[1, -1, 0], [-1, 2, -1], [0, -1, 1]])
        True
 	>>> print parallel.procID > 0 or numerix.allequal(solver.RHSvector, [0, 0, 0,])
        True
 	>>> eq = TransientTerm(coeff=1.) == DiffusionTerm(coeff=1., var=B) + 10. 
 	Traceback (most recent call last): 
 	    ... 
 	Exception: Terms with explicit Variables cannot mix with Terms with implicit Variables
 	>>> eq = DiffusionTerm(coeff=1., var=B) + 10. == 0 
 	>>> print eq 
 	(DiffusionTerm(coeff=[1.0], var=B) + 10.0)
 	>>> print eq._getVars()
 	[B]
 	>>> print (eq.term, eq.other)
        (DiffusionTerm(coeff=[1.0], var=B), 10.0)
 	>>> solver = eq._prepareLinearSystem(var=B, solver=None, boundaryConditions=(), dt=1.)
        >>> numpyMatrix = solver.matrix.getNumpyArray()
        >>> print parallel.procID > 0 or numerix.allequal(numpyMatrix, [[-1, 1, 0], [1, -2, 1], [0, 1, -1]])
        True
 	>>> print parallel.procID > 0 or numerix.allequal(solver.RHSvector, [-10, -10, -10]) 
 	True

 	>>> eq.solve(var=B)

        >>> m = Grid1D(nx=2)
        >>> A = CellVariable(mesh=m, name='A')
        >>> B = CellVariable(mesh=m, name='B')
        >>> C = CellVariable(mesh=m, name='C')        
        >>> DiffusionTerm().solve()
        Traceback (most recent call last):
            ...
        Exception: The solution variable needs to be specified
        >>> DiffusionTerm().solve(A)
        >>> DiffusionTerm(var=A).solve(A)
        >>> (DiffusionTerm(var=A) + DiffusionTerm())
        Traceback (most recent call last):
            ...
        Exception: Terms with explicit Variables cannot mix with Terms with implicit Variables
        >>> (DiffusionTerm(var=A) + DiffusionTerm(var=B)).solve()
        Traceback (most recent call last):
            ...
        Exception: The solution variable needs to be specified
        >>> (DiffusionTerm(var=A) + DiffusionTerm(var=B)).solve(A)
        >>> DiffusionTerm() & DiffusionTerm()
        Traceback (most recent call last):
            ...
        Exception: Different number of solution variables and equations.
        >>> DiffusionTerm(var=A) & DiffusionTerm()
        Traceback (most recent call last):
            ...
        Exception: Terms with explicit Variables cannot mix with Terms with implicit Variables
        >>> A = CellVariable(mesh=m, name='A', value=1)
        >>> B = CellVariable(mesh=m, name='B')
        >>> C = CellVariable(mesh=m, name='C')        
        >>> eq = (DiffusionTerm(coeff=1., var=A)) & (DiffusionTerm(coeff=2., var=B))
        >>> eq.cacheMatrix()
        >>> eq.cacheRHSvector()
        >>> eq.solve(solver=DefaultSolver())
        >>> numpyMatrix = eq.getMatrix().getNumpyArray()
        >>> print parallel.procID > 0 or numerix.allequal(numpyMatrix, [[-1, 1, 0, 0], [1, -1, 0, 0], [0, 0, -2, 2], [0, 0, 2, -2]])
        True
        >>> print eq.getRHSvector().getGlobalValue()
        [ 0.  0.  0.  0.]
        >>> print eq._getVars()
        [A, B]
        >>> DiffusionTerm(var=A) & DiffusionTerm(var=A)
        Traceback (most recent call last):
            ...
        Exception: Different number of solution variables and equations.
        >>> DiffusionTerm() & DiffusionTerm()
        Traceback (most recent call last):
            ...
        Exception: Different number of solution variables and equations.
        >>> (DiffusionTerm(var=A) & DiffusionTerm(var=B)).solve(A)
        Traceback (most recent call last):
            ...
        Exception: The solution variable should not be specified.
        >>> DiffusionTerm(var=A) & DiffusionTerm(var=B) & DiffusionTerm(var=B)
        Traceback (most recent call last):
            ...
        Exception: Different number of solution variables and equations.
        >>> (DiffusionTerm(var=A) & DiffusionTerm(var=B) & DiffusionTerm(var=C)).solve(solver=DefaultSolver())
        >>> (DiffusionTerm(var=A) & DiffusionTerm(var=B) & DiffusionTerm(var=C)).solve(A)
        Traceback (most recent call last):
            ...
        Exception: The solution variable should not be specified.
        >>> (DiffusionTerm(var=A) & (DiffusionTerm(var=B) + DiffusionTerm(var=C))).solve(A)
        Traceback (most recent call last):
            ...
        Exception: The solution variable should not be specified.
        >>> eq = (DiffusionTerm(coeff=1., var=A) + DiffusionTerm(coeff=2., var=B)) & (DiffusionTerm(coeff=2., var=B) + DiffusionTerm(coeff=3., var=C)) & (DiffusionTerm(coeff=3., var=C) + DiffusionTerm(coeff=1., var=A))
        >>> eq.cacheMatrix()
        >>> eq.cacheRHSvector()
        >>> eq.solve(solver=DefaultSolver())
        >>> numpyMatrix = eq.getMatrix().getNumpyArray()
        >>> print parallel.procID > 0 or numerix.allequal(numpyMatrix, [[-1, 1, -2, 2, 0, 0],
        ...                                                             [1, -1, 2, -2, 0, 0],
        ...                                                             [0, 0, -2, 2, -3, 3],
        ...                                                             [0, 0, 2, -2, 3, -3],
        ...                                                             [-1, 1, 0, 0, -3, 3],                
        ...                                                             [1, -1, 0, 0, 3, -3]])
        True
        >>> print eq.getRHSvector().getGlobalValue()
        [ 0.  0.  0.  0.  0.  0.]
        >>> print eq._getVars()
        [A, B, C]
            
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
