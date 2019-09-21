"""

This example solves the steady-state convection-diffusion equation as
described in :mod:`examples.diffusion.convection.exponential1D.mesh1D` but
uses the :class:`~fipy.terms.powerLawConvectionTerm.PowerLawConvectionTerm` rather than the
:class:`~fipy.terms.exponentialConvectionTerm.ExponentialConvectionTerm`.

>>> from fipy import CellVariable, Grid1D, DiffusionTerm, PowerLawConvectionTerm, Viewer
>>> from fipy.tools import numerix

>>> L = 10.
>>> nx = 1000
>>> mesh = Grid1D(dx = L / nx, nx = nx)

>>> valueLeft = 0.
>>> valueRight = 1.

>>> var = CellVariable(name = "concentration",
...                    mesh = mesh,
...                    value = valueLeft)

>>> var.constrain(valueLeft, mesh.facesLeft)
>>> var.constrain(valueRight, mesh.facesRight)

>>> diffCoeff = 1.
>>> convCoeff = (10.,)

>>> eq = (DiffusionTerm(coeff=diffCoeff)
...       + PowerLawConvectionTerm(coeff=convCoeff))

>>> eq.solve(var = var)

We test the solution against the analytical result:

>>> axis = 0
>>> x = mesh.cellCenters[axis]
>>> CC = 1. - numerix.exp(-convCoeff[axis] * x / diffCoeff)
>>> DD = 1. - numerix.exp(-convCoeff[axis] * L / diffCoeff)
>>> analyticalArray = CC / DD
>>> print(var.allclose(analyticalArray, rtol = 1e-2, atol = 1e-2))
1

If the problem is run interactively, we can view the result:

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

