#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "cellGradVariable.py"
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
 
__all__ = []

from fipy.variables.faceVariable import FaceVariable
from fipy.tools import numerix

class _FaceGradContributions(FaceVariable):
    """
    Test case

    >>> from fipy import *
    >>> m = Grid2D(nx=3, ny=3)
    >>> x, y = m.cellCenters
    >>> v = CellVariable(mesh=m, elementshape=(3,))
    >>> v[0] = x
    >>> v[1] = y
    >>> v[2] = x**2
    >>> out = _FaceGradContributions(v)
    >>> v0 = CellVariable(mesh=m, value=x)
    >>> v1 = CellVariable(mesh=m, value=y)
    >>> v2 = CellVariable(mesh=m, value=x**2)
    >>> print _FaceGradContributions(v).globalValue.shape
    (2, 3, 24)
    >>> print _FaceGradContributions(v0).globalValue.shape
    (2, 24)
    >>> print _FaceGradContributions(v0)
    [[ 0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.  -0.5  1.   2.
       2.5 -0.5  1.   2.   2.5 -0.5  1.   2.   2.5]
     [-0.5 -1.5 -2.5  0.5  1.5  2.5  0.5  1.5  2.5  0.5  1.5  2.5  0.   0.   0.
       0.   0.   0.   0.   0.   0.   0.   0.   0. ]]
    >>> print (_FaceGradContributions(v0).globalValue == out.globalValue[:,0]).all()
    True
    >>> print (_FaceGradContributions(v1).globalValue == out.globalValue[:,1]).all()
    True
    >>> print (_FaceGradContributions(v2).globalValue == out.globalValue[:,2]).all()
    True
    
    """

    def __init__(self, var):
        FaceVariable.__init__(self, mesh=var.mesh, elementshape=(var.mesh.dim,) + var.shape[:-1])
        self.var = self._requires(var)

    def _calcValue(self):
        faceValue = self.var.arithmeticFaceValue.numericValue
        return self.mesh._areaProjections[(slice(0,None,None),) + (numerix.newaxis,) * (len(faceValue.shape) - 1) + (slice(0,None,None),)] * faceValue[numerix.newaxis]
    
def _test(): 
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()
    
if __name__ == "__main__": 
    _test() 
