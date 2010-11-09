#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "coupledEquation.py"
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

from fipy.terms.term import Term
from fipy.matrices.pysparseMatrix import _CoupledPysparseMeshMatrix
from fipy.variables.cellVariable import CellVariable
from fipy.variables.coupledCellVariable import _CoupledCellVariable

class _CoupledEquation(Term):
    def __init__(self):
        Term.__init__(self)
        self.eqs = []
        
    def __repr__(self):
        return "\n  &\n".join((repr(eq) for eq in self.eqs))
        
    def _prepareLinearSystem(self, var, solver, boundaryConditions, dt):
        if var is not None:
            raise Exception("_CoupledEquation knows more about variables than you do")
            
        if boundaryConditions is not ():
            raise Exception("_CoupledEquation knows more about boundary conditions than you do")
            
        from fipy.solvers.pysparse import LinearLUSolver
        solver = LinearLUSolver()
#         solver = self.getDefaultSolver(solver)
        
        vars = set()
        for eq in self.eqs:
            vars = vars.union(set(eq.vars))
            
        if len(vars) != len(self.eqs):
            raise Exception("The number of equations is different from the number of unknowns")
        
        if os.environ.has_key('FIPY_DISPLAY_MATRIX'):
            if not hasattr(self, "_viewer"):
                from fipy.viewers.matplotlibViewer.matplotlibSparseMatrixViewer import MatplotlibSparseMatrixViewer
                Term._viewer = MatplotlibSparseMatrixViewer()

        bigVar = _CoupledCellVariable(vars=vars)
        bigRHSvector = []
        bigMatrix = []
        
        for eq in self.eqs:
            RHSvector = CellVariable(mesh=bigVar.getMesh())
            matrices = []
            for var in vars:
                solver_ij = eq._prepareLinearSystem(var=var, solver=None, boundaryConditions=(), dt=dt)
#                 print repr(var), eq
#                 print solver_ij.matrix
#                 print solver_ij.RHSvector
                matrices.append(solver_ij.matrix)
                RHSvector += solver_ij.RHSvector
            bigRHSvector.append(RHSvector)
            bigMatrix.append(matrices)
            
        bigMatrix = _CoupledPysparseMeshMatrix(mesh=bigVar.getMesh(), matrices=bigMatrix)
        bigRHSvector = _CoupledCellVariable(vars=bigRHSvector)
                
        solver._storeMatrix(var=bigVar, matrix=bigMatrix, RHSvector=bigRHSvector)
        
        if os.environ.has_key('FIPY_DISPLAY_MATRIX'):
            self._viewer.title = "%s %s" % (repr(bigVar), self.__class__.__name__)
            self._viewer.plot(matrix=bigMatrix, RHSvector=bigRHSvector)
            from fipy import raw_input
            raw_input()

        return solver
        
    def copy(self):
        eq = _CoupledEquation()
        eq.eqs.extend(self.eqs)
        
        return copy

    def __and__(self, other):
        dup = self.copy()
        dup &= other
        
        return dup
        
    def __iand__(self, other):
        self.eqs.append(other)
        
        return self

