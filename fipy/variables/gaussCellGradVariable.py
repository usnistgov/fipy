#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "gaussCellGradVariable.py"
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
 #  See the file "license.terms" for information on usage and  redistribution
 #  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #  
 # ###################################################################
 ##
 
__docformat__ = 'restructuredtext'

__all__ = []

from fipy.variables.cellVariable import CellVariable
from fipy.tools import numerix
from fipy.tools import inline
from fipy.variables.faceGradContributionsVariable import _FaceGradContributions

class _GaussCellGradVariable(CellVariable):
    """
    Test case for a vector cell variable

    >>> from fipy import *
    >>> m = Grid2D(nx=3, ny=3)
    >>> x, y = m.cellCenters
    >>> v = CellVariable(mesh=m, elementshape=(3,))
    >>> v[0] = x
    >>> v[1] = y
    >>> v[2] = x**2
    >>> v0 = CellVariable(mesh=m, value=x)
    >>> v1 = CellVariable(mesh=m, value=y)
    >>> v2 = CellVariable(mesh=m, value=x**2)
    >>> v.grad.globalValue.shape
    (2, 3, 9)
    >>> print v0.grad
    [[ 0.5  1.   0.5  0.5  1.   0.5  0.5  1.   0.5]
     [ 0.   0.   0.   0.   0.   0.   0.   0.   0. ]]
    >>> print (v0.grad.globalValue == v.grad.globalValue[:,0]).all()
    True
    >>> print (v1.grad.globalValue == v.grad.globalValue[:,1]).all()
    True
    >>> print (v2.grad.globalValue == v.grad.globalValue[:,2]).all()
    True
        
    """
    
    def __init__(self, var, name=''):
        CellVariable.__init__(self, mesh=var.mesh, name=name, elementshape=(var.mesh.dim,) + var.shape[:-1])
        self.var = self._requires(var)
        self.faceGradientContributions = _FaceGradContributions(self.var)


    def _calcValueInline(self, N, M, ids, orientations, volumes):
        val = self._array.copy()

        inline._runIterateElementInline("""
            ITEM(val, i, vec) = 0.;

            int k;
            for (k = 0; k < M; k++) {
                int id = ITEM(ids, i, &k);
                ITEM(val, i, vec) += ITEM(orientations, i, &k) * ITEM(areaProj, id, vec) * ITEM(faceValues, id, NULL);
            }

            ITEM(val, i, vec) /= ITEM(volumes, i, NULL);
        """,val = val,
            ids = numerix.array(numerix.MA.filled(ids, 0)),
            orientations = numerix.array(numerix.MA.filled(orientations, 0)),
            volumes = numerix.array(volumes),
            areaProj = numerix.array(self.mesh._areaProjections),
            faceValues = numerix.array(self.var.arithmeticFaceValue),
            M = M,
            ni = N, 
            shape=numerix.array(numerix.shape(val)))

        return self._makeValue(value = val)

    def _calcValueNoInline(self, N, M, ids, orientations, volumes):
        contributions = numerix.take(self.faceGradientContributions, ids, axis=-1)
        grad = numerix.array(numerix.sum(orientations * contributions, -2))
        return grad / volumes

    def _calcValue(self):
        if inline.doInline and self.var.rank == 0:
            return self._calcValueInline(N=self.mesh.numberOfCells, 
                                         M=self.mesh._maxFacesPerCell, 
                                         ids=self.mesh.cellFaceIDs, 
                                         orientations=self.mesh._cellToFaceOrientations, 
                                         volumes=self.mesh.cellVolumes)
        else:
            return self._calcValueNoInline(N=self.mesh.numberOfCells, 
                                           M=self.mesh._maxFacesPerCell, 
                                           ids=self.mesh.cellFaceIDs, 
                                           orientations=self.mesh._cellToFaceOrientations, 
                                           volumes=self.mesh.cellVolumes)
        

def _test(): 
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()
    
if __name__ == "__main__": 
    _test() 
