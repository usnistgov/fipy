"""

This input file solves a steady-state 1D diffusion problem as in
`./examples/diffusion/mesh1D.py`. The difference being that it uses a tensor for
the diffusion coefficient, even though the coefficient is isotropic.

>>> from fipy import Grid2D, CellVariable, DiffusionTerm, Viewer

>>> Lx = 20
>>> mesh = Grid2D(nx=20, ny=20)
>>> x, y = mesh.cellCenters

>>> valueLeft = 0.
>>> valueRight = 1.

>>> var = CellVariable(name = "solution variable",
...                    mesh = mesh,
...                    value = valueLeft)

>>> var.constrain(valueLeft, mesh.facesLeft)
>>> var.constrain(valueRight, mesh.facesRight)

>>> DiffusionTerm(coeff=(((1., 0.),
...                       (0., 1.)),)).solve(var)

>>> if __name__ == '__main__':
...     viewer = Viewer(vars=var).plot()

>>> analyticalArray = valueLeft + (valueRight - valueLeft) * x / Lx
>>> print(var.allclose(analyticalArray, atol = 0.025))
1
"""
from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy import input

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())

    input('finished')

