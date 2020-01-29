r"""

One can then solve the same problem as in
`examples/diffusion/steadyState/mesh1D/input.py` but with a periodic
mesh and no boundary conditions. The periodic mesh is used to simulate
periodic boundary conditions.

>>> from fipy import PeriodicGrid1D, CellVariable, TransientTerm, DiffusionTerm, Viewer

>>> nx = 50
>>> dx = 1.
>>> mesh = PeriodicGrid1D(nx = nx, dx = dx)

The variable is initially a line varying form `valueLeft` to `valueRight`.

>>> valueLeft = 0
>>> valueRight = 1
>>> x = mesh.cellCenters[0]

>>> Lx = nx * dx
>>> initialArray = valueLeft + (valueRight - valueLeft) * x / Lx
>>> var = CellVariable(name = "solution variable", mesh = mesh,
...                                                value = initialArray)

>>> from fipy import input
>>> if __name__ == '__main__':
...     viewer = Viewer(vars=var, datamin=0., datamax=1.)
...     viewer.plot()
...     input("press key to continue")


A `TransientTerm` is used to provide some fixed point, otherwise the
solver has no fixed value and can become unstable.

>>> eq = TransientTerm(coeff=1e-8) - DiffusionTerm()
>>> eq.solve(var=var, dt=1.)

>>> if __name__ == '__main__':
...     viewer.plot()

The result of the calculation will be the average value over the domain.

>>> print(var.allclose((valueLeft + valueRight) / 2., rtol = 1e-5))
1

"""
from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy import input

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())
    input("finished")
