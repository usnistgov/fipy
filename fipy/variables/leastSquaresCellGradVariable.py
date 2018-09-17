#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "leastSquaresCellGradVariable.py"
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

from fipy.variables.cellVariable import CellVariable
from fipy.tools import numerix

class _LeastSquaresCellGradVariable(CellVariable):
    """
    Look at CellVariable.leastSquarseGrad for documentation
     """
    def __init__(self, var, name = ''):
        CellVariable.__init__(self, mesh=var.mesh, name=name, rank=var.rank + 1)
        self.var = self._requires(var)

    @property
    def _neighborValue(self):
        return numerix.take(numerix.array(self.var), self.mesh._cellToCellIDs)

    def _calcValue(self):
        cellToCellDistances = self.mesh._cellToCellDistances
        cellNormals = self.mesh._cellNormals
        neighborValue = self._neighborValue
        value = numerix.array(self.var)
        cellDistanceNormals = cellToCellDistances * cellNormals

        N = self.mesh.numberOfCells
        M = self.mesh._maxFacesPerCell
        D = self.mesh.dim

        mat = numerix.zeros((D, D, M, N), 'd')

        ## good god! numpy.outer should have an axis argument!!!
        for i in range(D):
            for j in range(D):
                mat[i,j] = cellDistanceNormals[i] * cellDistanceNormals[j]

        mat = numerix.sum(mat, axis=2)

        vec = numerix.array(numerix.sum((neighborValue - value) * cellDistanceNormals, axis=1))

        if D == 1:
            vec[0] = vec[0] / mat[0, 0]
        elif D == 2:
            divisor = mat[0,0] * mat[1,1] - mat[0,1] * mat[1,0]
            gradx = (vec[0] * mat[1,1] - vec[1] * mat[1,0]) / divisor
            grady = (vec[1] * mat[0,0] - vec[0] * mat[0,1]) / divisor
            vec[0] = gradx
            vec[1] = grady
        else:
            ## very stuppy! numerix.linalg.solve should have an axis argument!!!
            for i in range(N):
                vec[...,i] = numerix.linalg.solve(mat[...,i],vec[...,i])

        return vec
