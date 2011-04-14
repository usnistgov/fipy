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

    def _verifyVar(self, var):
        if var is None:
            if len(self._vars) == 0:
                raise SolutionVariableRequiredError
            elif len(self._vars) == 1:
                return _BaseBinaryTerm._verifyVar(self, self._vars[0])
            else:
                return _BaseBinaryTerm._verifyVar(self, _CoupledCellVariable(self._vars))
        else:
            return var

    def _buildAndAddMatrices(self, var, SparseMatrix, boundaryConditions=(), dt=1.0, transientGeomCoeff=None, diffusionGeomCoeff=None):
        """Build matrices of constituent Terms and collect them
        
        Only called at top-level by `_prepareLinearSystem()`
        """
        if isinstance(var, _CoupledCellVariable) and len(var.vars) > 1:
            RHSvector = 0
            for v in var.vars:
                RHSvector += self._buildExplicitIfOtherVar(var=v, 
                                                           SparseMatrix=SparseMatrix, 
                                                           boundaryConditions=boundaryConditions, 
                                                           dt=dt,
                                                           transientGeomCoeff=transientGeomCoeff,
                                                           diffusionGeomCoeff=diffusionGeomCoeff)
                      
            matrix = SparseMatrix(mesh=var.mesh, numberOfVariables=len(var.vars), numberOfEquations=1)
        else:
            var, matrix, RHSvector = self._buildMatrix(var=var, 
                                                       SparseMatrix=SparseMatrix, 
                                                       boundaryConditions=boundaryConditions, 
                                                       dt=dt,
                                                       transientGeomCoeff=transientGeomCoeff,
                                                       diffusionGeomCoeff=diffusionGeomCoeff)
                                                       
            RHSvector += self._buildExplicitIfOtherVar(var=var, 
                                                       SparseMatrix=SparseMatrix, 
                                                       boundaryConditions=boundaryConditions, 
                                                       dt=dt,
                                                       transientGeomCoeff=transientGeomCoeff,
                                                       diffusionGeomCoeff=diffusionGeomCoeff)
                                                       
        return (var, matrix, RHSvector)

    def _buildMatrix(self, var, SparseMatrix,  boundaryConditions=(), dt=1.0, transientGeomCoeff=None, diffusionGeomCoeff=None):

        matrix = 0
        RHSvector = 0

        for term in (self.term, self.other):

            termVar, termMatrix, termRHSvector = term._buildMatrix(var,
                                                                   SparseMatrix,
                                                                   boundaryConditions=boundaryConditions,
                                                                   dt=dt,
                                                                   transientGeomCoeff=transientGeomCoeff,
                                                                   diffusionGeomCoeff=diffusionGeomCoeff)

            matrix += termMatrix
            RHSvector += termRHSvector
            
            term._buildCache(termMatrix, termRHSvector)

        if (os.environ.has_key('FIPY_DISPLAY_MATRIX')
            and os.environ['FIPY_DISPLAY_MATRIX'].lower() == "terms"): 
            self._viewer.title = "%s %s" % (var.name, repr(self))
            self._viewer.plot(matrix=matrix, RHSvector=RHSvector) 
            raw_input()

	return (var, matrix, RHSvector)

    def _buildExplicitIfOtherVar(self, var, SparseMatrix, boundaryConditions=(), dt=1.0, transientGeomCoeff=None, diffusionGeomCoeff=None):
        """return the residual vector if `var` is not `Term.var`, otherwise return 0
        """
        return (self.term._buildExplicitIfOtherVar(var=var, 
                                                   SparseMatrix=SparseMatrix, 
                                                   boundaryConditions=boundaryConditions, 
                                                   dt=dt,
                                                   transientGeomCoeff=transientGeomCoeff,
                                                   diffusionGeomCoeff=diffusionGeomCoeff) 
                + self.other._buildExplicitIfOtherVar(var=var, 
                                                      SparseMatrix=SparseMatrix, 
                                                      boundaryConditions=boundaryConditions, 
                                                      dt=dt,
                                                      transientGeomCoeff=transientGeomCoeff,
                                                      diffusionGeomCoeff=diffusionGeomCoeff))

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
