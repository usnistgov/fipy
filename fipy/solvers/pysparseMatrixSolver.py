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
 # of Standards and Technology, an agency of the Federal Government.
 # Pursuant to title 17 section 105 of the United States Code,
 # United States Code this software is not subject to copyright
 # protection, and this software is considered to be in the public domain.
 # FiPy is an experimental system.
 # NIST assumes no responsibility whatsoever for its use by whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 #
 # To the extent that NIST may hold copyright in countries other than the
 # United States, you are hereby granted the non-exclusive irrevocable and
 # unconditional right to print, publish, prepare derivative works and
 # distribute this software, in any medium, or authorize others to do so on
 # your behalf, on a royalty-free basis throughout the world.
 #
 # You may improve, modify, and create derivative works of the software or
 # any portion of the software, and you may copy and distribute such
 # modifications or works.  Modified works should carry a notice stating
 # that you changed the software and should note the date and nature of any
 # such change.  Please explicitly acknowledge the National Institute of
 # Standards and Technology as the original source.
 #
 # This software can be redistributed and/or modified freely provided that
 # any derivative works bear some notice that they are derived from it, and
 # any modified versions bear some notice that they have been modified.
 # ========================================================================
 #
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

__all__ = []

from fipy.solvers.solver import Solver
from fipy.matrices.pysparseMatrix import _PysparseMeshMatrix

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
