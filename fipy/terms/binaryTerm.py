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

from fipy.terms.baseBinaryTerm import _BaseBinaryTerm
from fipy.terms.explicitSourceTerm import _ExplicitSourceTerm
from fipy.variables.coupledCellVariable import _CoupledCellVariable
from fipy.terms import SolutionVariableRequiredError

class _BinaryTerm(_BaseBinaryTerm):

    @property
    def _buildExplcitIfOther(self):
        return True

    def _buildAndAddMatrices(self, var, SparseMatrix,  boundaryConditions=(), dt=1.0, transientGeomCoeff=None, diffusionGeomCoeff=None, buildExplicitIfOther=True):
        """Build matrices of constituent Terms and collect them

        Only called at top-level by `_prepareLinearSystem()`
        
        """

        matrix = SparseMatrix(mesh=var.mesh)
        RHSvector = 0

        for term in (self.term, self.other):
            
            tmpVar, tmpMatrix, tmpRHSvector = term._buildAndAddMatrices(var,
                                                                        SparseMatrix,
                                                                        boundaryConditions=boundaryConditions,
                                                                        dt=dt,
                                                                        transientGeomCoeff=transientGeomCoeff,
                                                                        diffusionGeomCoeff=diffusionGeomCoeff,
                                                                        buildExplicitIfOther=buildExplicitIfOther)

            matrix += tmpMatrix
            RHSvector += tmpRHSvector

            term._buildCache(tmpMatrix, tmpRHSvector)

        if (os.environ.has_key('FIPY_DISPLAY_MATRIX')
             and os.environ['FIPY_DISPLAY_MATRIX'].lower() == "terms"): 
             self._viewer.title = "%s %s" % (var.name, repr(self))
             self._viewer.plot(matrix=matrix, RHSvector=RHSvector) 
             raw_input()
             
        return (var, matrix, RHSvector)
    
    def _getDefaultSolver(self, solver, *args, **kwargs):
        for term in (self.term, self.other):
            defaultSolver = term._getDefaultSolver(solver, *args, **kwargs)
            if defaultSolver is not None:
                solver = defaultSolver
                
        return solver
        
    def __repr__(self):
        return '(' + repr(self.term) + ' + ' + repr(self.other) + ')'

    def __mul__(self, other):
        return other * self.term + other * self.other

    @property
    def _uncoupledTerms(self):
        return [self]

    def _getTransientGeomCoeff(self, var):
        return self._addNone(self.term._getTransientGeomCoeff(var), self.other._getTransientGeomCoeff(var))

    def _getDiffusionGeomCoeff(self, var):
        return self._addNone(self.term._getDiffusionGeomCoeff(var), self.other._getDiffusionGeomCoeff(var)) 

    __rmul__ = __mul__

    def _test(self):
        """
        >>> from fipy import *
        >>> m = Grid1D(nx=3)
        >>> v0 = CellVariable(mesh=m, value=0.)
        >>> v1 = CellVariable(mesh=m, value=1.)
        >>> eq = TransientTerm(var=v0) - DiffusionTerm(coeff=1., var=v0) - DiffusionTerm(coeff=2., var=v1)
        >>> var, matrix, RHSvector = eq._buildAndAddMatrices(var=v0, SparseMatrix=DefaultSolver()._matrixClass)
        >>> print var
        [ 0.  0.  0.]
        >>> print CellVariable(mesh=m, value=RHSvector).globalValue
        [ 0.  0.  0.]
        >>> print numerix.allequal(matrix.numpyArray, [[ 2, -1,  0],
        ...                                            [-1,  3, -1],
        ...                                            [ 0, -1,  2]])
        True
        >>> var, matrix, RHSvector = eq._buildAndAddMatrices(var=v1, SparseMatrix=DefaultSolver()._matrixClass)
        >>> print var
        [ 1.  1.  1.]
        >>> print CellVariable(mesh=m, value=RHSvector).globalValue
        [ 0.  0.  0.]
        >>> print numerix.allequal(matrix.numpyArray, [[ 2, -2,  0],
        ...                                            [-2,  4, -2],
        ...                                            [ 0, -2,  2]])
        True
        >>> print CellVariable(mesh=m, value=eq.justResidualVector(dt=1.)).globalValue
        [ 0.  0.  0.]
        
        >>> m = Grid1D(nx=6)
        >>> v0 = CellVariable(mesh=m, value=1.)
        >>> v1 = CellVariable(mesh=m, value=0.)
        >>> eq = TransientTerm(var=v0) - DiffusionTerm(coeff=1., var=v0) - DiffusionTerm(coeff=2., var=v1)
        >>> var, matrix, RHSvector = eq._buildAndAddMatrices(var=v0, SparseMatrix=DefaultSolver()._matrixClass) 
        >>> print var
        [ 1.  1.  1.  1.  1.  1.]
        >>> print CellVariable(mesh=m, value=RHSvector).globalValue
        [ 1.  1.  1.  1.  1.  1.]
        >>> print numerix.allequal(matrix.numpyArray, [[ 2,-1, 0, 0, 0, 0.],
        ...                                            [-1, 3,-1, 0, 0, 0.],
        ...                                            [ 0,-1, 3,-1, 0, 0.],
        ...                                            [ 0, 0,-1, 3,-1, 0.],
        ...                                            [ 0, 0, 0,-1, 3,-1.],
        ...                                            [ 0, 0, 0, 0,-1, 2.]])
        True
        >>> var, matrix, RHSvector = eq._buildAndAddMatrices(var=v1, SparseMatrix=DefaultSolver()._matrixClass) 
        >>> print var
        [ 0.  0.  0.  0.  0.  0.]
        >>> print CellVariable(mesh=m, value=RHSvector).globalValue
        [ 0.  0.  0.  0.  0.  0.]
        >>> print numerix.allequal(matrix.numpyArray, [[ 2,-2, 0, 0, 0, 0.],
        ...                                            [-2, 4,-2, 0, 0, 0.],
        ...                                            [ 0,-2, 4,-2, 0, 0.],
        ...                                            [ 0, 0,-2, 4,-2, 0.],
        ...                                            [ 0, 0, 0,-2, 4,-2.],
        ...                                            [ 0, 0, 0, 0,-2, 2.]])
        True
        >>> print CellVariable(mesh=m, value=eq.justResidualVector(dt=1.)).globalValue
        [ 0.  0.  0.  0.  0.  0.]

        >>> m = Grid1D(nx=3)
        >>> v0 = CellVariable(mesh=m, value=(0., 1., 2.))
        >>> v1 = CellVariable(mesh=m, value=(3., 4., 5.))
        >>> diffTerm = DiffusionTerm(coeff=1., var=v0)
        >>> eq00 = TransientTerm(var=v0) - diffTerm
        >>> eq0 = eq00 - DiffusionTerm(coeff=2., var=v1)
        >>> eq0.cacheMatrix()
        >>> diffTerm.cacheMatrix()
        >>> print CellVariable(mesh=m, value=eq0.justResidualVector(dt=1.)).globalValue
        [-3.  0.  3.]
        >>> eq0.solve(var=v0, solver=DummySolver())
        >>> print numerix.allequal(eq0.matrix.numpyArray, [[ 2, -1,  0],
        ...                                                [-1,  3, -1],
        ...                                                [ 0, -1,  2]])
        True
        >>> eq0.solve(var=v1, solver=DummySolver())
        >>> print numerix.allequal(eq0.matrix.numpyArray, [[ 2, -2,  0],
        ...                                                [-2,  4, -2],
        ...                                                [ 0, -2,  2]])
        True
        >>> ## This currectly returns None because we lost the handle to the DiffusionTerm when it's negated.
        >>> print diffTerm.matrix 
        None

        Testing solution for one variable in a multi-variable equation.

        >>> from fipy import *
        >>> L = 1.
        >>> nx = 3
        >>> m = Grid1D(nx=nx, dx=L / nx)
        >>> x = m.cellCenters[0]
        >>> v0 = CellVariable(mesh=m, value=0., name='v0')
        >>> v0.constrain(0., where=m.facesLeft)
        >>> v0.constrain(L, where=m.facesRight)
        >>> v1 = CellVariable(mesh=m, value=-x**2, name='v1')
        >>> v1.constrain(0.,  where=m.facesLeft)
        >>> v1.constrain(-L,  where=m.facesRight)
        >>> (DiffusionTerm(var=v0) + DiffusionTerm(var=v1)).solve(v0)
        >>> print numerix.allclose(v0, -v1)
        True
        
        """


def _test(): 
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
