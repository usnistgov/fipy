#!/usr/bin/env python

##
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "gmshinput.py"
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
`./examples/diffusion/steadyState/mesh1D/input.py`. In order to test
the non-orthogonality error, this uses a SkewedGrid2D, which is a
Grid2D with each interior vertex moved in a random direction.

"""

if __name__ == '__main__':
    import sys
    import os

    from fipy import SkewedGrid2D, CellVariable, Viewer, DiffusionTerm
    from fipy.tools import numerix

    valueLeft = 0.
    valueRight = 1.

    mesh = SkewedGrid2D(dx = 1.0, dy = 1.0, nx = 20, ny = 20, rand = 0.1)

    var = CellVariable(name = "solution variable",
                       mesh = mesh,
                       value = valueLeft)

    viewer = Viewer(vars = var)

    var.constrain(valueLeft, mesh.facesLeft)
    var.constrain(valueRight, mesh.facesRight)

    DiffusionTerm().solve(var)

    varArray = numerix.array(var)
    x = mesh.cellCenters[0]
    analyticalArray = valueLeft + (valueRight - valueLeft) * x / 20
    errorArray = varArray - analyticalArray
    errorVar = CellVariable(name = 'absolute error',
                            mesh = mesh,
                            value = abs(errorArray))
    errorViewer = Viewer(vars = errorVar)

    NonOrthoVar = CellVariable(name = "non-orthogonality",
                               mesh = mesh,
                               value = mesh._nonOrthogonality)
    NOViewer = Viewer(vars = NonOrthoVar)
    viewer.plot()
    NOViewer.plot()

    raw_input("finished")
