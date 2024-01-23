# -*- coding: utf-8 -*-


r"""Solve the steady-state convection-diffusion equation in one dimension.

This example solves the steady-state convection-diffusion equation
given by

.. math::

   \nabla \cdot \left(D \nabla \phi + \vec{u} \phi \right) = 0

with coefficients :math:`D = 1` and :math:`\vec{u} = 10\hat{\imath}`, or

>>> diffCoeff = 1.
>>> convCoeff = (10.,)

We define a 1D mesh

.. index::
   single: Grid1D

>>> from fipy import CellVariable, Grid1D, DiffusionTerm, ExponentialConvectionTerm, Viewer
>>> from fipy.tools import numerix

>>> L = 10.
>>> nx = 10
>>> mesh = Grid1D(dx=L / nx, nx=nx)

>>> valueLeft = 0.
>>> valueRight = 1.

The solution variable is initialized to ``valueLeft``:

>>> var = CellVariable(mesh=mesh, name="variable")

and impose the boundary conditions

.. math::

   \phi = \begin{cases}
   0& \text{at $x = 0$,} \\
   1& \text{at $x = L$,}
   \end{cases}

with

>>> var.constrain(valueLeft, mesh.facesLeft)
>>> var.constrain(valueRight, mesh.facesRight)

The equation is created with the :class:`~fipy.terms.diffusionTerm.DiffusionTerm` and
:class:`~fipy.terms.exponentialConvectionTerm.ExponentialConvectionTerm`. The scheme used by the convection term
needs to calculate a Péclet number and thus the diffusion term
instance must be passed to the convection term.

>>> eq = (DiffusionTerm(coeff=diffCoeff)
...       + ExponentialConvectionTerm(coeff=convCoeff))

More details of the benefits and drawbacks of each type of convection
term can be found in :ref:`sec:NumericalSchemes`.
Essentially, the :class:`~fipy.terms.exponentialConvectionTerm.ExponentialConvectionTerm` and :class:`~fipy.terms.powerLawConvectionTerm.PowerLawConvectionTerm` will
both handle most types of convection-diffusion cases, with the
:class:`~fipy.terms.powerLawConvectionTerm.PowerLawConvectionTerm` being more efficient.

We solve the equation

>>> eq.solve(var=var)

and test the solution against the analytical result

.. math::

   \phi = \frac{1 - \exp(-u_x x / D)}{1 - \exp(-u_x L / D)}

or

.. index::
   single: exp

>>> axis = 0
>>> x = mesh.cellCenters[axis]
>>> CC = 1. - numerix.exp(-convCoeff[axis] * x / diffCoeff)
>>> DD = 1. - numerix.exp(-convCoeff[axis] * L / diffCoeff)
>>> analyticalArray = CC / DD
>>> print(var.allclose(analyticalArray))
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

