#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "leastSquaresCellGradVariable.py"
 #                                    created: 12/18/03 {2:28:00 PM} 
 #                                last update: 1/7/09 {9:24:05 AM} 
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
 #  See the file "license.terms" for information on usage and  redistribution
 #  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #  
 # ###################################################################
 ##
 
from fipy.variables.vectorCellVariable import VectorCellVariable
from fipy.tools import numerix

class _LeastSquaresCellGradVariable(VectorCellVariable):
    """
    Look at CellVariable.getLeastSquarseGrad() for documentation
    """
    def __init__(self, var, name = ''):
        VectorCellVariable.__init__(self, mesh = var.getMesh(), name = name)
        self.var = self._requires(var)

    def _getNeighborValue(self, ):
        return numerix.take(numerix.array(self.var), self.mesh._getCellToCellIDs())

    def _calcValue(self):
        cellToCellDistances = self.mesh._getCellToCellDistances()
        cellNormals = self.mesh._getCellNormals()
        neighborValue = self._getNeighborValue()
        value = numerix.array(self.var)
        cellDistanceNormals = cellToCellDistances[...,numerix.newaxis] * cellNormals

        N = self.mesh.getNumberOfCells()
        M = self.mesh._getMaxFacesPerCell()
        D = self.mesh.getDim()

        mat = numerix.zeros((N, M, D, D), 'd')

        ## good god! numpy.outer should have and axis argument!!!
        for i in range(D):
            for j in range(D):
                mat[:,:,i,j] = cellDistanceNormals[:, :, i] * cellDistanceNormals[:, :, j]

        mat = numerix.sum(mat, axis=1)

        # (neighborValue - value[...,numerix.newaxis])[...,numerix.newaxis] * cellDistanceNormals
        # produces erroneous results (mask and sometimes value) with NumPy 1.2 and 1.3.
        # Bug http://scipy.org/scipy/numpy/ticket/979 has been filed.

        vec = numerix.sum((neighborValue - value[...,numerix.newaxis])[...,numerix.newaxis] * cellDistanceNormals, axis=1)

        if D == 1:
            vec[:,0] = vec[:,0] / mat[:, 0, 0]
        elif D == 2:
            divisor = mat[:,0,0] * mat[:,1,1] - mat[:,0,1] * mat[:,1,0]
            gradx = (vec[:,0] * mat[:,1,1] - vec[:,1] * mat[:,0,1]) / divisor
            grady = (vec[:,1] * mat[:,0,0] - vec[:,0] * mat[:,1,0]) / divisor
            vec[:,0] = gradx
            vec[:,1] = grady
        else:
            ## very stuppy! numerix.linalg.solve should have an axis argument!!!        
            for i in range(N):
                vec[i] = numerix.linalg.solve(mat[i],vec[i])

        return vec

