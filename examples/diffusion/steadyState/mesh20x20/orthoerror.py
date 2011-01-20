#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "orthoerror.py"
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

"""

This test file generates lots of different SkewedGrid2D meshes, each with a different non-orthogonality,
and runs a 1D diffusion problem on them all. It ocmputes the RMS non-orthogonality and the RMS error
for each mesh and displays them in a graph, allowing the relationship of error to non-orthogonality to be investigated.
"""

if __name__ == '__main__':
    
    import sys
    import os

    from fipy import *

    valueLeft = 0.
    valueRight = 1.

    meshList = []
    RMSNonOrthoList = []
    RMSErrorList = []

    for i in range(1, 501):
        meshList = meshList + [SkewedGrid2D(dx = 1.0, dy = 1.0, nx = 20, ny = 20, rand = (0.001 * i))]

    for mesh in meshList:
        var = CellVariable(name = "solution variable",
                           mesh = mesh,
                           value = valueLeft)

        var.constrain(valueLeft, mesh.facesLeft)
        var.constrain(valueRight, mesh.facesRight)

        DiffusionTerm().solve(var)

        varArray = array(var)
        x = mesh.cellCenters[0]
        analyticalArray = valueLeft + (valueRight - valueLeft) * x / 20
        errorArray = varArray - analyticalArray
        nonOrthoArray = mesh._nonOrthogonality
        RMSError = (add.reduce(errorArray * errorArray) / len(errorArray)) ** 0.5
        RMSNonOrtho = (add.reduce(nonOrthoArray * nonOrthoArray) / len(nonOrthoArray)) ** 0.5

        RMSNonOrthoList += [RMSNonOrtho]
        RMSErrorList += [RMSError]

    import pylab
    pylab.plot(RMSNonOrthoList, RMSErrorList, 'ro')
    pylab.show()



