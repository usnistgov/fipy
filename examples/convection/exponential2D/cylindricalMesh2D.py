r"""

This example solves the steady-state cylindrical convection-diffusion equation
given by:

.. math::

   \nabla \cdot \left(D \nabla \phi + \vec{u} \phi \right) = 0

with coefficients :math:`D = 1` and :math:`\vec{u} = (10,)`, or

>>> diffCoeff = 1.
>>> convCoeff = ((10.,), (0.,))

We define a 2D cylindrical mesh representing an annulus. The mesh is a
pseudo 1D mesh, but is a good test case for the
:func:`~fipy.meshes.factoryMeshes.CylindricalGrid2D` mesh.

.. index::
   single: Grid1D

>>> from fipy import CellVariable, CylindricalGrid2D, DiffusionTerm, ExponentialConvectionTerm, Viewer
>>> from fipy.tools import numerix

>>> r0 = 1.
>>> r1 = 2.
>>> nr = 100
>>> mesh = CylindricalGrid2D(dr=(r1 - r0) / nr, dz=1., nr=nr, nz=1) + ((r0,),)

The solution variable is initialized to ``valueLeft``:

>>> valueLeft = 0.
>>> valueRight = 1.
>>> var = CellVariable(mesh=mesh, name = "variable")

and impose the boundary conditions

.. math::

   \phi = \begin{cases}
   0& \text{at $r = r_0$,} \\
   1& \text{at $r = r_1$,}
   \end{cases}

with

>>> var.constrain(valueLeft, mesh.facesLeft)
>>> var.constrain(valueRight, mesh.facesRight)

The equation is created with the :class:`~fipy.terms.diffusionTerm.DiffusionTerm` and
:class:`~fipy.terms.exponentialConvectionTerm.ExponentialConvectionTerm`.

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

   \phi = \exp{\frac{u}{D} \left(r_1 - r\right)} \left( \frac{ \Ei{\frac{u r_0}{D}} - \Ei{\frac{u r}{D}} }{ \Ei{\frac{u r_0}{D}} - \Ei{\frac{u r_1}{D}} } \right)

or

.. index::
   single: exp

>>> axis = 0
>>> try:
...     from scipy.special import expi # doctest: +SCIPY
...     U = convCoeff[0][0]
...     r = mesh.cellCenters[axis]
...     AA = numerix.exp(U / diffCoeff * (r1 - r))
...     BB = expi(U * r0 / diffCoeff) - expi(U * r / diffCoeff) # doctest: +SCIPY
...     CC = expi(U * r0 / diffCoeff) - expi(U * r1 / diffCoeff) # doctest: +SCIPY
...     analyticalArray = AA * BB / CC # doctest: +SCIPY
... except ImportError:
...     print("The SciPy library is unavailable. It is required for testing purposes.")

>>> print(var.allclose(analyticalArray, atol=1e-3)) # doctest: +SCIPY
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


