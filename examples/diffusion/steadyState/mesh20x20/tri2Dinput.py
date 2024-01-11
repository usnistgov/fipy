"""

This input file again solves a 2D diffusion problem on a triangular mesh.
We increase the solver tolerance from the default :math:`10^{-5}` in order
to achieve a good solution.

    >>> eq = DiffusionTerm()
    >>> solver = eq.getDefaultSolver(tolerance=1e-10)
    >>> eq.solve(var, solver=solver)

The result is again tested in the same way:

    >>> Lx = nx * dx
    >>> x = mesh.cellCenters[0]
    >>> analyticalArray = valueLeft + (valueRight - valueLeft) * x / Lx
    >>> print(var.allclose(analyticalArray, rtol = 1e-8))
    True

"""
from __future__ import unicode_literals

from fipy import input
from fipy import CellVariable, Tri2D, DiffusionTerm, Viewer

nx = 20
ny = 20

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

