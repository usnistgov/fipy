##from fipy.tools.profiler.profiler import Profiler
##from fipy.tools.profiler.profiler import calibrate_profiler

"""

This input file again solves a 1D diffusion problem as in
`./examples/diffusion/steadyState/mesh1D/input.py`. The difference
being that the mesh size is given by

    >>> nx = 50
    >>> ny = 50

The result is again tested in the same way:

    >>> DiffusionTerm().solve(var)
    >>> Lx = nx * dx
    >>> x = mesh.cellCenters[0]
    >>> analyticalArray = valueLeft + (valueRight - valueLeft) * x / Lx
    >>> print(var.allclose(analyticalArray, atol = 1e-7))
    1

"""
from __future__ import unicode_literals

from fipy import input
from fipy import CellVariable, Tri2D, DiffusionTerm, Viewer

nx = 50
ny = 50

dx = 1.

valueLeft = 0.
valueRight = 1.

mesh = Tri2D(dx = dx, nx = nx, ny = ny)

var = CellVariable(name = "solution variable",
                   mesh = mesh,
                   value = valueLeft)

var.constrain(valueLeft, mesh.facesLeft)
var.constrain(valueRight, mesh.facesRight)

if __name__ == '__main__':
    DiffusionTerm().solve(var)
    viewer = Viewer(vars = var)
    viewer.plot()
    input("finished")

