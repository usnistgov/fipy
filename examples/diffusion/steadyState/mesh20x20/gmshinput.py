##from fipy.tools.profiler.profiler import Profiler
##from fipy.tools.profiler.profiler import calibrate_profiler

"""

This input file again solves a 1D diffusion problem as in
`./examples/diffusion/steadyState/mesh1D/input.py`. In order to test
the non-orthogonality error, this uses a `SkewedGrid2D`, which is a
`Grid2D` with each interior vertex moved in a random direction.

"""
from __future__ import division
from __future__ import unicode_literals

from fipy import input
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

    input("finished")
