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

from fipy.terms.baseBinaryTerm import _BaseBinaryTerm
from fipy.variables.coupledCellVariable import _CoupledCellVariable
from fipy.variables.cellVariable import CellVariable
from fipy.tools import numerix
from fipy.terms import SolutionVariableNumberError
from fipy.terms import AlternativeMethodInBaseClass

class _CoupledBinaryTerm(_BaseBinaryTerm):
    """
    Test to ensure that _getTransientGeomCoeff and _getDiffusionGeomCoeff return sensible results for coupled equations.

    >>> from fipy import *
    >>> m = Grid1D(nx=1)
    >>> v0 = CellVariable(mesh=m)
    >>> v1 = CellVariable(mesh=m)
    >>> eq0 = TransientTerm(1, var=v0) == DiffusionTerm(2, var=v1)
    >>> eq1 = TransientTerm(3, var=v1) == DiffusionTerm(4, var=v0)
    >>> eq = eq0 & eq1
    >>> print eq._getTransientGeomCoeff(v0)
    None
    >>> print eq._getDiffusionGeomCoeff(v1)
    None
    >>> tranCoeff = eq._getUncoupledTerms()[0]._getTransientGeomCoeff(v0)
    >>> print parallel.procID > 0 or numerix.allequal(tranCoeff, [1])
    True
    >>> diffCoeff = eq._getUncoupledTerms()[1]._getDiffusionGeomCoeff(v0)
    >>> print parallel.procID > 0 or numerix.allequal(diffCoeff, [[-8, -8]])
    True
    
    """
    def __init__(self, term, other):
        _BaseBinaryTerm.__init__(self, term, other)
        if len(self._getVars()) < len(self._getUncoupledTerms()):
            raise SolutionVariableNumberError

    def _getUncoupledTerms(self):
        return self.term._getUncoupledTerms() + self.other._getUncoupledTerms()

    def _verifyVar(self, var):
        if var is not None:
            raise Exception, 'The solution variable should not be specified.'

        if len(self._getVars()) != len(self._getUncoupledTerms()):
            raise SolutionVariableNumberError

        return _BaseBinaryTerm._verifyVar(self, _CoupledCellVariable(self._getVars()))
    
    def _buildMatrix(self, var, SparseMatrix,  boundaryConditions=(), dt=1.0, transientGeomCoeff=None, diffusionGeomCoeff=None):
        """
        Offset tests

        >>> from fipy import *
        >>> m = Grid1D(nx=3)
        >>> v0 = CellVariable(mesh=m, value=0.)
        >>> v1 = CellVariable(mesh=m, value=1.)        
        >>> eq0 = TransientTerm(var=v0) - DiffusionTerm(coeff=1., var=v0) - DiffusionTerm(coeff=2., var=v1)
        >>> eq1 = TransientTerm(var=v1) - DiffusionTerm(coeff=3., var=v0) - DiffusionTerm(coeff=4., var=v1) 
        >>> eq = eq0 & eq1
        >>> var = eq._verifyVar(None)
        >>> solver = DefaultSolver()
        >>> var, matrix, RHSvector = eq._buildMatrix(var=var, SparseMatrix=DefaultSolver()._matrixClass) 
        >>> print var.globalValue
        [ 0.  0.  0.  1.  1.  1.]
        >>> print RHSvector.globalValue
        [ 0.  0.  0.  1.  1.  1.]
        >>> print numerix.allequal(matrix.asTrilinosMeshMatrix().numpyArray,
        ...                        [[2, -1, 0, 2, -2, 0],
        ...                         [-1, 3, -1, -2, 4, -2],
        ...                         [0, -1, 2, 0, -2, 2],
        ...                         [3, -3, 0, 5, -4, 0],
        ...                         [-3, 6, -3, -4, 9, -4],                
        ...                         [0, -3, 3, 0, -4, 5]])
        True

        >>> m = Grid1D(nx=6)
        >>> v0 = CellVariable(mesh=m, value=0.)
        >>> v1 = CellVariable(mesh=m, value=1.)        
        >>> eq0 = TransientTerm(var=v0) - DiffusionTerm(coeff=1., var=v0) - DiffusionTerm(coeff=2., var=v1)
        >>> eq1 = TransientTerm(var=v1) - DiffusionTerm(coeff=3., var=v0) - DiffusionTerm(coeff=4., var=v1) 
        >>> eq = eq0 & eq1
        >>> var = eq._verifyVar(None)
        >>> solver = DefaultSolver()
        >>> var, matrix, RHSvector = eq._buildMatrix(var=var, SparseMatrix=DefaultSolver()._matrixClass) 
        >>> print var.globalValue
        [ 0.  0.  0.  0.  0.  0.  1.  1.  1.  1.  1.  1.]
        >>> print RHSvector.globalValue
        [ 0.  0.  0.  0.  0.  0.  1.  1.  1.  1.  1.  1.]
        >>> print numerix.allequal(matrix.asTrilinosMeshMatrix().numpyArray,
        ...                        [[ 2, -1,  0,  0,  0,  0,  2, -2,  0,  0,  0,  0],
        ...                         [-1,  3, -1,  0,  0,  0, -2,  4, -2,  0,  0,  0],
        ...                         [ 0, -1,  3, -1,  0,  0,  0, -2,  4, -2,  0,  0],
        ...                         [ 0,  0, -1,  3, -1,  0,  0,  0, -2,  4, -2,  0],
        ...                         [ 0,  0,  0, -1,  3, -1,  0,  0,  0, -2,  4, -2],
        ...                         [ 0,  0,  0,  0, -1,  2,  0,  0,  0,  0, -2,  2],
        ...                         [ 3, -3,  0,  0,  0,  0,  5, -4,  0,  0,  0,  0],
        ...                         [-3,  6, -3,  0,  0,  0, -4,  9, -4,  0,  0,  0],
        ...                         [ 0, -3,  6, -3,  0,  0,  0, -4,  9, -4,  0,  0],
        ...                         [ 0,  0, -3,  6, -3,  0,  0,  0, -4,  9, -4,  0],
        ...                         [ 0,  0,  0, -3,  6, -3,  0,  0,  0, -4,  9, -4],
        ...                         [ 0,  0,  0,  0, -3,  3,  0,  0,  0,  0, -4,  5]])
        True
        
        >>> m = Grid1D(nx=3)
        >>> v0 = CellVariable(mesh=m, value=0.)
        >>> v1 = CellVariable(mesh=m, value=1.)
        >>> diffTerm = DiffusionTerm(coeff=1., var=v0)
        >>> eq00 = TransientTerm(var=v0) - DiffusionTerm(coeff=1., var=v0)
        >>> eq0 = eq00 - DiffusionTerm(coeff=2., var=v1)
        >>> eq1 = TransientTerm(var=v1) - DiffusionTerm(coeff=3., var=v0) - DiffusionTerm(coeff=4., var=v1) 
        >>> eq0.cacheMatrix()
        >>> diffTerm.cacheMatrix()
        >>> (eq0 & eq1).solve()
        >>> print numerix.allequal(eq0.getMatrix().asTrilinosMeshMatrix().numpyArray,
        ...                        [[ 0,  0,  0,  2, -2,  0],
        ...                         [ 0,  0,  0, -2,  4, -2],
        ...                         [ 0,  0,  0,  0, -2,  2],
        ...                         [ 0,  0,  0,  0,  0,  0],
        ...                         [ 0,  0,  0,  0,  0,  0],
        ...                         [ 0,  0,  0,  0,  0,  0]])
        True
        >>> ## This currectly returns None because we lost the handle to the DiffusionTerm when it's negated.
        >>> print diffTerm.getMatrix() 
        None
        
        """

        numberOfCells = var.mesh.numberOfCells
        numberOfVariables = len(self._getVars())
        
        matrix = 0
        RHSvectorsJ = []

        for i, uncoupledTerm in enumerate(self._getUncoupledTerms()):

            RHSvector = 0

            for j, tmpVar in enumerate(self._getVars()):

                class OffsetSparseMatrix(SparseMatrix):
                    def __init__(self, mesh, bandwidth=0, sizeHint=None, numberOfVariables=numberOfVariables):
                        SparseMatrix.__init__(self, mesh=mesh, bandwidth=bandwidth, sizeHint=sizeHint, numberOfVariables=numberOfVariables)

                    def put(self, vector, id1, id2):
                        SparseMatrix.put(self, vector, id1 + numberOfCells * i, id2 + numberOfCells * j)

                    def addAt(self, vector, id1, id2):
                        SparseMatrix.addAt(self, vector, id1 + numberOfCells * i, id2 + numberOfCells * j)

                    def addAtDiagonal(self, vector):
                        if type(vector) in [type(1), type(1.)]:
                            tmp = numerix.zeros((numberOfCells,), 'd')
                            tmp[:] = vector
                            SparseMatrix.addAtDiagonal(self, tmp)
                        else:
                            SparseMatrix.addAtDiagonal(self, vector)

                tmpVar, tmpMatrix, tmpRHSvector = uncoupledTerm._buildMatrix(tmpVar,
                                                                             OffsetSparseMatrix,
                                                                             boundaryConditions=(),
                                                                             dt=dt,
                                                                             transientGeomCoeff=uncoupledTerm._getTransientGeomCoeff(tmpVar),
                                                                             diffusionGeomCoeff=uncoupledTerm._getDiffusionGeomCoeff(tmpVar))

                uncoupledTerm._buildCache(tmpMatrix, tmpRHSvector)

                RHSvector += tmpRHSvector
                matrix += tmpMatrix

            RHSvectorsJ += [CellVariable(value=RHSvector, mesh=var.mesh)]

        RHSvector = _CoupledCellVariable(RHSvectorsJ)

	return (var, matrix, RHSvector)

    def __repr__(self):
        return '(' + repr(self.term) + ' & ' + repr(self.other) + ')'

    def _getDefaultSolver(self, solver, *args, **kwargs):
        if _BaseBinaryTerm._getDefaultSolver(self, solver, *args, **kwargs) is not None:
            raise AlternativeMethodInBaseClass('getDefaultSolver()')

        if solver and not solver._canSolveAsymmetric():
            import warnings
            warnings.warn("%s cannot solve assymetric matrices" % solver)
        from fipy.solvers import DefaultAsymmetricSolver
        return solver or DefaultAsymmetricSolver(*args, **kwargs)    

    def _calcVars(self):
        """
        This method returns the equations variables ordered by, transient terms,
        diffusion terms and other terms. Currently, this won't cure all of
        coupled equations ills, but it fixes the majority of issues, mainly with
        preconditioning. Some tests to make sure it works correctly.

        >>> from fipy import *
        >>> m = Grid1D(nx=2)
        >>> v0 = CellVariable(mesh=m, name='v0')
        >>> v1 = CellVariable(mesh=m, name='v1')
        >>> v2 = CellVariable(mesh=m, name='v2')
        >>> ConvectionTerm = CentralDifferenceConvectionTerm
        >>> eq0 = TransientTerm(var=v0) + VanLeerConvectionTerm(var=v0) + DiffusionTerm(var=v1)
        >>> eq1 = TransientTerm(var=v2) + ConvectionTerm(var=v2) + DiffusionTerm(var=v2) + ConvectionTerm(var=v1) + ImplicitSourceTerm(var=v1)
        >>> eq2 = ImplicitSourceTerm(var=v1) + 1 + ImplicitSourceTerm(var=v0) + 1 + DiffusionTerm(var=v1)
        >>> print (eq0 & eq1 & eq2)._getVars()
        [v0, v2, v1]
        >>> print (eq0 & eq2 & eq1)._getVars()
        [v0, v1, v2]
        >>> eq0 =  DiffusionTerm(var=v1) + TransientTerm(var=v0) + VanLeerConvectionTerm(var=v0)
        >>> print (eq0 & eq2 & eq1)._getVars()
        [v0, v1, v2]
        >>> print (eq2 & eq0 & eq1)._getVars()
        [v1, v0, v2]
        >>> print (eq2 & eq0 & eq1)([v1, v2, v0])._getVars()
        [v1, v2, v0]
        >>> print (eq2 & eq0 & eq1)([v1, v2, v0, v2])._getVars()
  	Traceback (most recent call last): 
 	    ... 
 	SolutionVariableNumberError: Different number of solution variables and equations.
        >>> print (eq2 & eq0 & eq1)([v1, v2, 1])._getVars()
  	Traceback (most recent call last): 
 	    ... 
 	Exception: Variable not in previously defined variables for this coupled equation.
        >>> print (eq2 & eq0 & eq1)([v1, v2, v1])._getVars()
 	Traceback (most recent call last): 
 	    ... 
 	SolutionVariableNumberError: Different number of solution variables and equations.
        >>> print (eq2 & eq0 & eq1)([v1, v2])._getVars()
 	Traceback (most recent call last): 
 	    ... 
 	SolutionVariableNumberError: Different number of solution variables and equations.

        """
    
        ## set() is used to force comparison by reference rather than value
        unorderedVars = _BaseBinaryTerm._calcVars(self)
        uncoupledTerms = self._getUncoupledTerms()
        
        if len(unorderedVars) == len(uncoupledTerms):
            unorderedVars = set(unorderedVars)
            orderedVars = [None] * len(uncoupledTerms)

            for fnc in (lambda index, term: term._getTransientVars(),
                        lambda index, term: term._getDiffusionVars(),
                        lambda index, term: list(unorderedVars)):
                for index, term in enumerate(uncoupledTerms):
                    if orderedVars[index] is None:
                        _vars = fnc(index, term)
                        if  _vars != [] and _vars[0] in unorderedVars:
                            orderedVars[index] = _vars[0]
                            unorderedVars.remove(_vars[0])
            
            return orderedVars
        else:
            ## Constituent _CoupledBinaryTerms don't necessarily have the same
            ## number of equations and variables so ordering is unnecessary.
            return unorderedVars

    def __call__(self, _vars):
        _vars = list(_vars)

        if len(_vars) != len(self._getVars()) or len(set(_vars)) != len(self._getVars()):
            raise SolutionVariableNumberError

        for var in _vars:
            if var not in set(self._getVars()):
                raise Exception, 'Variable not in previously defined variables for this coupled equation.'

        self._vars = _vars
        return self

def _test(): 
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()

