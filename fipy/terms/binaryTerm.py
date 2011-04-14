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

class _BinaryTerm(_BaseBinaryTerm):

    def _verifyVars(self, var):
        if var is None:
            if len(self._vars) == 0:
                raise SolutionVariableRequiredError

        return var, self._vars
        ##     elif len(self._vars) == 1:
        ##         return _BaseBinaryTerm._verifyVar(self, self._vars[0])
        ##     else:
        ##         return _BaseBinaryTerm._verifyVar(self, _CoupledCellVariable(self._vars))
        ## else:
        ##     return var

    def _getMatrixClass(self, solver):
        return solver._matrixClass
##      from fipy.matrices.offsetSparseMatrix import OffsetSparseMatrix
##      return OffsetSparseMatrix(SparseMatrix=solver._matrixClass,
##                                numberOfVariables=len(self._vars) or 1,
##                                                                                              
        
    def _buildAndAddMatrices(self, solutionVar, equationVars, SparseMatrix, boundaryConditions=(), dt=1.0, transientGeomCoeff=None, diffusionGeomCoeff=None):
        """Build matrices of constituent Terms and collect them
        
        Only called at top-level by `_prepareLinearSystem()`

        We want _BinaryTerm to return a rectangular matrix.

        >>> from fipy import *
        >>> m = Grid1D(nx=3)
        >>> v0 = CellVariable(mesh=m, value=0.)
        >>> v1 = CellVariable(mesh=m, value=1.)
        >>> eq = TransientTerm(var=v0) - DiffusionTerm(coeff=1., var=v0) - DiffusionTerm(coeff=2., var=v1)
        >>> var, matrix, RHSvector = eq._buildAndAddMatrices(var=eq._verifyVar(None), SparseMatrix=eq._getMatrixClass(DefaultSolver()))
        >>> print var.globalValue
        [ 1.  1.  1.  0.  0.  0.]
        >>> print CellVariable(mesh=m, value=RHSvector).globalValue
        [ 0.  0.  0.]
        >>> print numerix.allequal(matrix.numpyArray,
        ...                        [[ 2, -2,  0,  2, -1,  0],
        ...                         [-2,  4, -2, -1,  3, -1],
        ...                         [ 0, -2,  2,  0, -1,  2]])
        True
        
        >>> m = Grid1D(nx=6)
        >>> v0 = CellVariable(mesh=m, value=1.)
        >>> v1 = CellVariable(mesh=m, value=0.)
        >>> eq = TransientTerm(var=v0) - DiffusionTerm(coeff=1., var=v0) - DiffusionTerm(coeff=2., var=v1)
        >>> var, matrix, RHSvector = eq._buildAndAddMatrices(var=eq._verifyVar(None), SparseMatrix=eq._getMatrixClass(DefaultSolver())) 
        >>> print var.globalValue
        [ 0.  0.  0.  0.  0.  0.  1.  1.  1.  1.  1.  1.]
        >>> print CellVariable(mesh=m, value=RHSvector).globalValue
        [ 1.  1.  1.  1.  1.  1.]
        >>> print numerix.allequal(matrix.numpyArray,
        ...                        [[ 2, -2,  0,  0,  0,  0,  2, -1,  0,  0,  0,  0],
        ...                         [-2,  4, -2,  0,  0,  0, -1,  3, -1,  0,  0,  0],
        ...                         [ 0, -2,  4, -2,  0,  0,  0, -1,  3, -1,  0,  0],
        ...                         [ 0,  0, -2,  4, -2,  0,  0,  0, -1,  3, -1,  0],
        ...                         [ 0,  0,  0, -2,  4, -2,  0,  0,  0, -1,  3, -1],
        ...                         [ 0,  0,  0,  0, -2,  2,  0,  0,  0,  0, -1,  2]])
        True
                
        >>> m = Grid1D(nx=3)
        >>> v0 = CellVariable(mesh=m, value=(0., 1., 2.))
        >>> v1 = CellVariable(mesh=m, value=(3., 4., 5.))
        >>> eq00 = TransientTerm(var=v0) - DiffusionTerm(coeff=1., var=v0)
        >>> eq0 = eq00 - DiffusionTerm(coeff=2., var=v1)
        >>> eq0.cacheMatrix()
        >>> diffTerm.cacheMatrix()
        >>> print eq0.justResidualVector()
        [-3.  0.  3.]
        >>> print numerix.allequal(eq0.matrix.numpyArray,
        ...                        [[ 2, -2,  0,  2, -1,  0],
        ...                         [-2,  4, -2, -1,  3, -1],
        ...                         [ 0, -2,  2,  0, -1,  2]])
        True
        >>> ## This currectly returns None because we lost the handle to the DiffusionTerm when it's negated.
        >>> print diffTerm.matrix 
        None
        
        """

        matrix = SparseMatrix(mesh=solutionVar.mesh)
        RHSvector = 0

        for term in (self.term, self.other):
        
            tmpVar, tmpMatrix, tmpRHSvector = term._buildAndAddMatrices(solutionVar,
                                                                        equationVars,
                                                                        SparseMatrix,
                                                                        boundaryConditions=boundaryConditions,
                                                                        dt=dt,
                                                                        transientGeomCoeff=transientGeomCoeff,
                                                                        diffusionGeomCoeff=diffusionGeomCoeff)

            matrix += tmpMatrix
            RHSvector += tmpRHSvector

            term._buildCache(tmpMatrix, tmpRHSvector)

        if (os.environ.has_key('FIPY_DISPLAY_MATRIX')
             and os.environ['FIPY_DISPLAY_MATRIX'].lower() == "terms"): 
             self._viewer.title = "%s %s" % (var.name, repr(self))
             self._viewer.plot(matrix=matrix, RHSvector=RHSvector) 
             raw_input()

        return (solutionVar, matrix, RHSvector)

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

def _test(): 
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
