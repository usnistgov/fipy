r"""Solve the steady-state convection-diffusion equation with a constant source.

Like :mod:`examples.convection.exponential1D.mesh1D`
this example solves a steady-state convection-diffusion equation, but adds a constant source,
:math:`S_0 = 1`, such that

.. math::

   \nabla \cdot \left(D \nabla \phi + \vec{u} \phi \right) + S_0 = 0.

>>> diffCoeff = 1.
>>> convCoeff = (10.,)
>>> sourceCoeff = 1.

We define a 1D mesh

.. index::
   single: Grid1D

>>> from fipy import CellVariable, Grid1D, DiffusionTerm, ExponentialConvectionTerm, DefaultAsymmetricSolver, Viewer
>>> from fipy.tools import numerix

>>> nx = 1000
>>> L = 10.
>>> mesh = Grid1D(dx=L / 1000, nx=nx)

>>> valueLeft = 0.
>>> valueRight = 1.

The solution variable is initialized to ``valueLeft``:

.. index::
   single: CellVariable

>>> var = CellVariable(name="variable", mesh=mesh)

and impose the boundary conditions

.. math::

   \phi = \begin{cases}
   0& \text{at $x = 0$,} \\
   1& \text{at $x = L$,}
   \end{cases}

with

>>> var.constrain(valueLeft, mesh.facesLeft)
>>> var.constrain(valueRight, mesh.facesRight)

We define the convection-diffusion equation with source

>>> eq = (DiffusionTerm(coeff=diffCoeff)
...       + ExponentialConvectionTerm(coeff=convCoeff)
...       + sourceCoeff)

.. index::
   single: DefaultAsymmetricSolver

>>> eq.solve(var=var,
...          solver=DefaultAsymmetricSolver(tolerance=1.e-15, iterations=10000))

and test the solution against the analytical result:

.. math::

   \phi = -\frac{S_0 x}{u_x}
   + \left(1 + \frac{S_0 x}{u_x}\right)\frac{1 - \exp(-u_x x / D)}{1 - \exp(-u_x L / D)}

or

.. index::
   single: exp

>>> axis = 0
>>> x = mesh.cellCenters[axis]
>>> AA = -sourceCoeff * x / convCoeff[axis]
>>> BB = 1. + sourceCoeff * L / convCoeff[axis]
>>> CC = 1. - numerix.exp(-convCoeff[axis] * x / diffCoeff)
>>> DD = 1. - numerix.exp(-convCoeff[axis] * L / diffCoeff)
>>> analyticalArray = AA + BB * CC / DD
>>> print(var.allclose(analyticalArray, rtol=1e-4, atol=1e-4))
1

If the problem is run interactively, we can view the result:

.. index::
   pair: module; fipy.viewers

>>> if __name__ == '__main__':
...     viewer = Viewer(vars=var)
...     viewer.plot()
"""
from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy import input

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())

    input('finished')

