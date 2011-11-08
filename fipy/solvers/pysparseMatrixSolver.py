#!/usr/bin/env python

## 
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "pysparseSolver.py"
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

__all__ = []

from fipy.solvers.solver import Solver
from fipy.matrices.pysparseMatrix import _PysparseMeshMatrix

from fipy.tools.decorators import getsetDeprecated
 
class _PysparseMatrixSolver(Solver):

    """
    A class consolidating methods for solver packages which use
    `_PysparseMeshMatrix` for their matrix class.

    Subclasses have a `_solve_` method, which is called by `_solve`. Typically,
    `_solve_` returns the new value of `self.var` to `_solve` and solve sets the
    var accordingly.

    A solution function `solveFnc`, usually of the form `solve(A, x, b)`, is
    implemented in most leaf-node child classes.

    .. attention:: This class is abstract. Always create one of its subclasses.
    """

    solveFnc = None

    @getsetDeprecated
    def _getMatrixClass(self):
        return self._matrixClass

    @property
    def _matrixClass(self):
        return _PysparseMeshMatrix
         
    def _solve(self):
        """
        Call `_solve_` for the new value of `self.var`.

        In certain cases, `_solve_` won't return anything, e.g. 
        `fipy.solvers.pysparse.linearLUSolver`. In these cases, we preserve the
        value of `self.var.numericValue`.
        """

        if self.var.mesh.communicator.Nproc > 1:
            raise Exception("%ss cannot be used with multiple processors" \
                            % self.__class__)
        
        array = self.var.numericValue
        newArr = self._solve_(self.matrix, array, self.RHSvector)

        if newArr is not None:
            array = newArr

        factor = self.var.unit.factor

        if factor != 1:
            array /= self.var.unit.factor

        self.var[:] = array  
