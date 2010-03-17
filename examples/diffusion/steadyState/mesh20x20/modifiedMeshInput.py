#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "ttri2Dinput.py"
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

##from fipy.tools.profiler.profiler import Profiler
##from fipy.tools.profiler.profiler import calibrate_profiler

"""

This input file again solves a 1D diffusion problem as in
`./examples/diffusion/steadyState/mesh1D/input.py`. The difference
being that it uses a triangular mesh loaded in using the GmshImporter.

The result is again tested in the same way:

    >>> DiffusionTerm().solve(var, boundaryConditions = boundaryConditions)
    >>> Lx = 20
    >>> x = mesh.getCellCenters()[0]
    >>> analyticalArray = valueLeft + (valueRight - valueLeft) * x / Lx
    >>> print var.allclose(analyticalArray, atol = 0.025)
    1

    >>> max(mesh._getNonOrthogonality()) < 0.51
    True

Note that this test case will only work if you run it by running the
main FiPy test suite. If you run it directly from the directory it is
in it will not be able to find the mesh file.

"""

import sys

from fipy import *

valueLeft = 0.
valueRight = 1.

import os.path
mesh = GmshImporter2D(os.path.join(os.path.split(__file__)[0], 'modifiedMesh.msh'))

##    "%s/%s" % (sys.__dict__['path'][0], "examples/diffusion/steadyState/mesh20x20/modifiedMesh.msh"))

var = CellVariable(name = "solution variable",
                   mesh = mesh,
                   value = valueLeft)

exteriorFaces = mesh.getExteriorFaces()
xFace = mesh.getFaceCenters()[0]
boundaryConditions = (FixedValue(exteriorFaces & (xFace ** 2 < 0.000000000000001), valueLeft),
                      FixedValue(exteriorFaces & ((xFace - 20) ** 2 < 0.000000000000001), valueRight))
                      

if __name__ == '__main__':
    DiffusionTerm().solve(var, boundaryConditions = boundaryConditions)
    viewer = Viewer(vars = var)
    viewer.plot()
    varArray = array(var)
    x = mesh.getCellCenters()[0]
    analyticalArray = valueLeft + (valueRight - valueLeft) * x / 20
    errorArray = varArray - analyticalArray
    errorVar = CellVariable(name = "absolute error",
                   mesh = mesh,
                   value = abs(errorArray))
    errorViewer = Viewer(vars = errorVar)
    errorViewer.plot()

    NonOrthoVar = CellVariable(name = "non-orthogonality",
                               mesh = mesh,
                               value = mesh._getNonOrthogonality())
    NOViewer = Viewer(vars = NonOrthoVar)


    NOViewer.plot()
    raw_input("finished")
