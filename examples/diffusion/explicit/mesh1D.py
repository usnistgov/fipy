r"""

This input file again solves a 1D diffusion problem as in
:mod:`examples.diffusion.steadyState.mesh1D`,
the difference being that this transient example is solved explicitly.

We create a 1D mesh:

>>> from fipy import CellVariable, Grid1D, TransientTerm, ExplicitDiffusionTerm, Viewer
>>> from fipy.tools import numerix

>>> nx = 100
>>> dx = 1.
>>> mesh = Grid1D(dx = dx, nx = nx)

and we initialize a :class:`~fipy.variables.cellVariable.CellVariable` to
``initialValue``:

>>> valueLeft = 0.
>>> initialValue = 1.
>>> var = CellVariable(
...     name = "concentration",
...     mesh = mesh,
...     value = initialValue)

The transient diffusion equation

.. math::

   \frac{\partial \phi}{\partial t} = \nabla \cdot (D \nabla \phi)

is represented by a :class:`~fipy.terms.transientTerm.TransientTerm` and an
:class:`~fipy.terms.explicitDiffusionTerm.ExplicitDiffusionTerm`.

We take the diffusion coefficient :math:`D = 1`

>>> diffusionCoeff = 1.

We build the equation:

>>> eq = TransientTerm() == ExplicitDiffusionTerm(coeff = diffusionCoeff)

and the boundary conditions:

>>> var.constrain(valueLeft, mesh.facesLeft)

In this case, many steps have to be taken to reach equilibrium.  A loop is
required to execute the necessary time steps:

>>> timeStepDuration = 0.1
>>> steps = 100
>>> from builtins import range
>>> for step in range(steps):
...     eq.solve(var=var, dt=timeStepDuration)

The analytical solution for this transient diffusion problem is given
by :math:`\phi = \erf(x/2\sqrt{D t})`.

The result is tested against the expected profile:

>>> Lx = nx * dx
>>> x = mesh.cellCenters[0]
>>> t = timeStepDuration * steps
>>> epsi = x / numerix.sqrt(t * diffusionCoeff)
>>> from scipy.special import erf # doctest: +SCIPY
>>> analyticalArray = erf(epsi/2) # doctest: +SCIPY
>>> print(var.allclose(analyticalArray, atol = 2e-3)) # doctest: +SCIPY
1

If the problem is run interactively, we can view the result:

>>> if __name__ == '__main__':
...     viewer = Viewer(vars = (var,))
...     viewer.plot()
"""
from __future__ import unicode_literals

__docformat__ = 'restructuredtext'

from fipy import input

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())
    input('finished')
