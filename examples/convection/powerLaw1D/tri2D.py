"""

This example solves the steady-state convection-diffusion equation as described in
:mod:`examples.convection.exponential1D.mesh1D` but uses the
:class:`~fipy.terms.powerLawConvectionTerm.PowerLawConvectionTerm` rather than the :class:`~fipy.terms.exponentialConvectionTerm.ExponentialConvectionTerm` instantiator.

>>> from fipy import CellVariable, Grid2D, DiffusionTerm, PowerLawConvectionTerm, DefaultAsymmetricSolver, Viewer
>>> from fipy.tools import numerix

>>> L = 10.
>>> nx = 1000
>>> mesh = Grid2D(dx = L / nx, nx = nx)

>>> valueLeft = 0.
>>> valueRight = 1.

>>> var = CellVariable(name = "concentration",
...                    mesh = mesh,
...                    value = valueLeft)

>>> var.constrain(valueLeft, mesh.facesLeft)
>>> var.constrain(valueRight, mesh.facesRight)

>>> diffCoeff = 1.
>>> convCoeff = (10., 0.)

>>> eq = (DiffusionTerm(coeff=diffCoeff)
...       + PowerLawConvectionTerm(coeff=convCoeff))

>>> eq.solve(var=var,
...          solver=DefaultAsymmetricSolver(tolerance=1.e-15, iterations=2000))

The analytical solution test for this problem is given by:

>>> axis = 0
>>> x = mesh.cellCenters[axis]
>>> CC = 1. - numerix.exp(-convCoeff[axis] * x / diffCoeff)
>>> DD = 1. - numerix.exp(-convCoeff[axis] * L / diffCoeff)
>>> analyticalArray = CC / DD
>>> print(var.allclose(analyticalArray, rtol = 1e-2, atol = 1e-2))
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


