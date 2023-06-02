"""

This input file again solves a 1D diffusion problem as in
:mod:`examples.diffusion.mesh1D`.
The difference being that the mesh is two dimensional.

The result is again tested in the same way:

    >>> DiffusionTerm().solve(var)
    >>> Lx = nx * dx
    >>> x = mesh.cellCenters[0]
    >>> analyticalArray = valueLeft + (valueRight - valueLeft) * x / Lx
    >>> print(var.allclose(analyticalArray, rtol = 1e-9))
    1

"""
from __future__ import unicode_literals

__docformat__ = 'restructuredtext'

from fipy import input
from fipy import CellVariable, Grid2D, DiffusionTerm, Viewer

nx = 50
ny = 50

dx = 1.

valueLeft = 0.
valueRight = 1.

mesh = Grid2D(dx = dx, nx = nx, ny = ny)

var = CellVariable(name = "solution variable",
                   mesh = mesh,
                   value = valueLeft)

var.constrain(valueLeft, mesh.facesLeft)
var.constrain(valueRight, mesh.facesRight)

if __name__ == '__main__':
    DiffusionTerm().solve(var)

    viewer = Viewer(vars=var, datamin=0., datamax=1.)
    viewer.plot()
    input("finished")

