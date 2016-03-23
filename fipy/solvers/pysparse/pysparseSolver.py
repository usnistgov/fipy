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

import os
from fipy.solvers.pysparseMatrixSolver import _PysparseMatrixSolver

__all__ = ["PysparseSolver"]

class PysparseSolver(_PysparseMatrixSolver):
    """
    The base `pysparseSolver` class.

    .. attention:: This class is abstract. Always create one of its subclasses.
    """
    def __init__(self, *args, **kwargs):
        if self.__class__ is PysparseSolver:
            raise NotImplementedError("can't instantiate abstract base class")

        super(PysparseSolver, self).__init__(*args, **kwargs)

    def _solve_(self, L, x, b):
        """
        `_solve_` is only for use by solvers which may use
        preconditioning. If you are writing a solver which
        doesn't use preconditioning, this must be overridden.

        :Parameters:
            - `L`: a `fipy.matrices.pysparseMatrix._PysparseMeshMatrix`.
            - `x`: a `numpy.ndarray`.
            - `b`: a `numpy.ndarray`.
        """

        A = L.matrix

        if self.preconditioner is None:
            P = None
        else:
            P, A = self.preconditioner._applyToMatrix(A)

        info, iter, relres = self.solveFnc(A, b, x, self.tolerance,
                                           self.iterations, P)

        self._raiseWarning(info, iter, relres)

        if 'FIPY_VERBOSE_SOLVER' in os.environ:
            from fipy.tools.debug import PRINT
            PRINT('iterations: %d / %d' % (iter, self.iterations))

            if info < 0:
                PRINT('failure', self._warningList[info].__class__.__name__)
            PRINT('relres:', relres)

    def _solve(self):

        if self.var.mesh.communicator.Nproc > 1:
            raise Exception("PySparse solvers cannot be used with multiple processors")

        array = self.var.numericValue.ravel()

        from fipy.terms import SolutionVariableNumberError

        if ((self.matrix == 0)
            or (self.matrix.matrix.shape[0] != self.matrix.matrix.shape[1])
            or (self.matrix.matrix.shape[0] != len(array))):

            raise SolutionVariableNumberError

        self._solve_(self.matrix, array, self.RHSvector)
        factor = self.var.unit.factor
        if factor != 1:
            array /= self.var.unit.factor

        self.var[:] = array.reshape(self.var.shape)
