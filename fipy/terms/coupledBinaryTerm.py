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
from fipy.tools import isin

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
    >>> print eq._getUncoupledTerms()[0]._getTransientGeomCoeff(v0)
    [ 1.]
    >>> print eq._getUncoupledTerms()[1]._getDiffusionGeomCoeff(v0)
    [[-8. -8.]]
    
    """
    def __init__(self, term, other):
        _BaseBinaryTerm.__init__(self, term, other)
        if len(self._getVars()) < len(self._getUncoupledTerms()):
            raise Exception, 'Different number of solution variables and equations.'

    def _getUncoupledTerms(self):
        return self.term._getUncoupledTerms() + self.other._getUncoupledTerms()

    def _verifyVar(self, var):
        if var is not None:
            raise Exception, 'The solution variable should not be specified.'

        if len(self._getVars()) != len(self._getUncoupledTerms()):
            raise Exception, 'Different number of solution variables and equations.'

        return _BaseBinaryTerm._verifyVar(self, _CoupledCellVariable(self.__getOrderedVars()))
    
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
        >>> var, matrix, RHSvector = eq._buildMatrix(var=var, SparseMatrix=DefaultSolver()._getMatrixClass()) 
        >>> print var.getGlobalValue()
        [ 0.  0.  0.  1.  1.  1.]
        >>> print RHSvector.getGlobalValue()
        [ 0.  0.  0.  1.  1.  1.]
        >>> print numerix.allequal(matrix.asTrilinosMeshMatrix().getNumpyArray(),
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
        >>> var, matrix, RHSvector = eq._buildMatrix(var=var, SparseMatrix=DefaultSolver()._getMatrixClass()) 
        >>> print var.getGlobalValue()
        [ 0.  0.  0.  0.  0.  0.  1.  1.  1.  1.  1.  1.]
        >>> print RHSvector.getGlobalValue()
        [ 0.  0.  0.  0.  0.  0.  1.  1.  1.  1.  1.  1.]
        >>> print numerix.allequal(matrix.asTrilinosMeshMatrix().getNumpyArray(),
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

        
        
        """

        numberOfCells = var.getMesh().getNumberOfCells()
        numberOfVariables = len(self._getVars())
        
        matrix = 0
        RHSvectorsJ = []

        for i, uncoupledTerm in enumerate(self._getUncoupledTerms()):

            RHSvector = 0

            for j, tmpVar in enumerate(self.__getOrderedVars()):

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

                RHSvector += tmpRHSvector
                matrix += tmpMatrix

            RHSvectorsJ += [CellVariable(value=RHSvector, mesh=var.getMesh())]

        RHSvector = _CoupledCellVariable(RHSvectorsJ)

	return (var, matrix, RHSvector)

    def __repr__(self):
        return '(' + repr(self.term) + ' & ' + repr(self.other) + ')'

    def _getDefaultSolver(self, solver, *args, **kwargs):
        if _BaseBinaryTerm._getDefaultSolver(self, solver, *args, **kwargs) is not None:
            raise AssertionError, 'An alternate _getDefaultSolver() is defined in a base class'

        if solver and not solver._canSolveAsymmetric():
            import warnings
            warnings.warn("%s cannot solve assymetric matrices" % solver)
        from fipy.solvers import DefaultAsymmetricSolver
        return solver or DefaultAsymmetricSolver(*args, **kwargs)    

    def _getOrderedVarsForTests(self):
        return self.__getOrderedVars()
        
    def __getOrderedVars(self):
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
        >>> print (eq0 & eq1 & eq2)._getOrderedVarsForTests()
        [v0, v2, v1]
        >>> print (eq0 & eq2 & eq1)._getOrderedVarsForTests()
        [v0, v1, v2]
        >>> eq0 =  DiffusionTerm(var=v1) + TransientTerm(var=v0) + VanLeerConvectionTerm(var=v0)
        >>> print (eq0 & eq2 & eq1)._getOrderedVarsForTests()
        [v0, v1, v2]
        >>> print (eq2 & eq0 & eq1)._getOrderedVarsForTests()
        [v1, v0, v2]
        >>> print (eq2 & eq0 & eq1)([v1, v2, v0])._getOrderedVarsForTests()
        [v1, v2, v0]
        >>> print (eq2 & eq0 & eq1)([v1, v2, v0, v2])._getOrderedVarsForTests()
  	Traceback (most recent call last): 
 	    ... 
 	Exception: Different number of solution variables and equations.
        >>> print (eq2 & eq0 & eq1)([v1, v2, 1])._getOrderedVarsForTests()
  	Traceback (most recent call last): 
 	    ... 
 	Exception: Variable not in previously defined variables for this coupled equation.
        >>> print (eq2 & eq0 & eq1)([v1, v2, v1])._getOrderedVarsForTests()
 	Traceback (most recent call last): 
 	    ... 
 	Exception: Different number of solution variables and equations.
        >>> print (eq2 & eq0 & eq1)([v1, v2])._getOrderedVarsForTests()
 	Traceback (most recent call last): 
 	    ... 
 	Exception: Different number of solution variables and equations.

        """
        if not hasattr(self, '_orderedVars'):
            uncoupledTerms = self._getUncoupledTerms()
            unorderedVars = self._getVars()[:]
            orderedVars = [None] * len(uncoupledTerms)


            for fnc in (lambda index, term: term._getTransientVars(),
                        lambda index, term: term._getDiffusionVars(),
                        lambda index, term: unorderedVars):

                for index, term in enumerate(uncoupledTerms):
                    if orderedVars[index] is None:
                        _vars = fnc(index, term)
                        if  _vars != [] and isin(_vars[0], unorderedVars):
                            orderedVars[index] = _vars[0]
                            unorderedVars = [u for u in unorderedVars if u is not _vars[0]]
                            ## unorderedVars.remove(_vars[0]) Does not compare references, but values. 
                            
            self._orderedVars = orderedVars
                
        return self._orderedVars

    def __call__(self, _orderedVars):
        _orderedVars = list(_orderedVars)

        if len(_orderedVars) != len(self._getVars()):
            raise Exception, 'Different number of solution variables and equations.'
        import fipy.tools

        _orderedVars = fipy.tools.uniqueList(_orderedVars)

        for var in _orderedVars:
            if not isin(var, self._getVars()):            
                raise Exception, 'Variable not in previously defined variables for this coupled equation.'

        if len(_orderedVars) != len(self._getVars()):
            raise Exception, 'Different number of solution variables and equations.'
        self._orderedVars = _orderedVars

        return self
        
               
def _test(): 
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()

