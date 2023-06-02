r"""

This example solves the steady-state convection-diffusion equation as described in
:mod:`examples.convection.exponential1D.mesh1D` but uses a
:class:`~fipy.meshes.tri2D.Tri2D` mesh.

Here the axes are reversed (``nx = 1``, ``ny = 1000``) and

.. math::

   \vec{u} = (0, 10)

>>> from fipy import CellVariable, Tri2D, DiffusionTerm, ExponentialConvectionTerm, DefaultAsymmetricSolver, Viewer
>>> from fipy.tools import numerix

>>> L = 10.
>>> nx = 1
>>> ny = 1000
>>> mesh = Tri2D(dx = L / ny, dy = L / ny, nx = nx, ny = ny)

>>> valueBottom = 0.
>>> valueTop = 1.

>>> var = CellVariable(name = "concentration",
...                    mesh = mesh,
...                    value = valueBottom)

>>> var.constrain(valueBottom, mesh.facesBottom)
>>> var.constrain(valueTop, mesh.facesTop)

>>> diffCoeff = 1.
>>> convCoeff = numerix.array(((0.,), (10.,)))

>>> eq = (DiffusionTerm(coeff=diffCoeff)
...       + ExponentialConvectionTerm(coeff=convCoeff))

>>> eq.solve(var = var,
...          solver=DefaultAsymmetricSolver(iterations=10000))

The analytical solution test for this problem is given by:

>>> axis = 1
>>> y = mesh.cellCenters[axis]
>>> CC = 1. - numerix.exp(-convCoeff[axis] * y / diffCoeff)
>>> DD = 1. - numerix.exp(-convCoeff[axis] * L / diffCoeff)
>>> analyticalArray = CC / DD

>>> print(var.allclose(analyticalArray, rtol = 1e-6, atol = 1e-6))
1

>>> if __name__ == '__main__':
...     viewer = Viewer(vars = var)
...     viewer.plot()
"""
from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy import input

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())

    input('finished')

