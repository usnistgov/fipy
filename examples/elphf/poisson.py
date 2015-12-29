#!/usr/bin/env python

##
 # ###################################################################
 #  FiPy - Python-based phase field solver
 #
 #  FILE: "input1DpoissonRightCharge.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this software is not subject to copyright
 # protection and is in the public domain.  FiPy is an experimental
 # system.  NIST assumes no responsibility whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 #
 # This software can be redistributed and/or modified freely
 # provided that any derivative works bear some notice that they are
 # derived from it, and any modified versions bear some notice that
 # they have been modified.
 # ========================================================================
 #
 # ###################################################################
 ##

r"""
A simple 1D example to test the setup of the Poisson equation.

>>> from fipy import CellVariable, Grid1D, DiffusionTerm, Viewer
>>> from fipy.tools import numerix

>>> nx = 200
>>> dx = 0.01
>>> L = nx * dx
>>> mesh = Grid1D(dx = dx, nx = nx)

The dimensionless Poisson equation is

.. math::

   \nabla\cdot\left(\epsilon\nabla\phi\right) = -\rho = -\sum_{j=1}^n z_j C_j

where :math:`\phi` is the electrostatic potential, :math:`\epsilon` is the
permittivity, :math:`\rho` is the charge density, :math:`C_j` is the
concentration of the :math:`j^\text{th}` component, and :math:`z_j` is the
valence of the :math:`j^\text{th}` component.

We will be solving for the electrostatic potential

>>> potential = CellVariable(mesh = mesh, name = 'phi', value = 0.)
>>> permittivity = 1.

We examine a fixed distribution of electrons with :math:`z_{\text{e}^{-}} = -1`.

>>> class ComponentVariable(CellVariable):
...     def __init__(self, mesh, value = 0., name = '',
...                  standardPotential = 0., barrier = 0.,
...                  diffusivity = None, valence = 0, equation = None):
...         CellVariable.__init__(self, mesh = mesh,
...                               value = value, name = name)
...         self.standardPotential = standardPotential
...         self.barrier = barrier
...         self.diffusivity = diffusivity
...         self.valence = valence
...         self.equation = equation
...
...     def copy(self):
...         return self.__class__(mesh = self.mesh,
...                               value = self.value,
...                               name = self.name,
...                               standardPotential =
...                                   self.standardPotential,
...                               barrier = self.barrier,
...                               diffusivity = self.diffusivity,
...                               valence = self.valence,
...                               equation = self.equation)

Since we're only interested in a single species, electrons, we could
simplify the following, but since we will in general be studying multiple
components, we explicitly allow for multiple substitutional species and
multiple interstitial species:

>>> interstitials = [
...     ComponentVariable(mesh = mesh, name = 'e-', valence = -1)]
>>> substitutionals = []

Because Poisson's equation admits an infinite number of potential profiles,
we must constrain the solution by fixing the potential at one point:

>>> potential.constrain(0., mesh.facesLeft)

>>> charge = 0.
>>> for Cj in interstitials + substitutionals:
...     charge += Cj * Cj.valence

>>> potential.equation = DiffusionTerm(coeff = permittivity) \
...                      + charge == 0

First, we obtain a uniform charge distribution by setting a uniform
concentration of electrons :math:`C_{\text{e}^{-}} = 1`.

>>> interstitials[0].setValue(1.)

and we solve for the electrostatic potential

>>> potential.equation.solve(var = potential)

This problem has the analytical solution

.. math::

   \psi(x) = \frac{x^2}{2} - 2x

We verify that the correct equilibrium is attained

>>> x = mesh.cellCenters[0]
>>> analyticalArray = (x**2)/2 - 2*x

>>> print potential.allclose(analyticalArray, rtol = 2e-5, atol = 2e-5)
1

If we are running the example interactively, we view the result

>>> if __name__ == '__main__':
...     viewer = Viewer(vars = (charge, potential))
...     viewer.plot()
...     raw_input("Press any key to continue...")

Next, we segregate all of the electrons to right side of the domain

.. math::

   C_{\text{e}^{-}} =
   \begin{cases}
       0& \text{for $x \le L/2$,} \\
       1& \text{for $x > L/2$.}
   \end{cases}

>>> x = mesh.cellCenters[0]
>>> interstitials[0].setValue(0.)
>>> interstitials[0].setValue(1., where=x > L / 2.)

and again solve for the electrostatic potential

>>> potential.equation.solve(var = potential)

which now has the analytical solution

.. math::

   \psi(x) =
   \begin{cases}
       -x& \text{for $x \le L/2$,} \\
       \frac{(x-1)^2}{2} - x& \text{for $x > L/2$.}
   \end{cases}

We verify that the correct equilibrium is attained

>>> analyticalArray = numerix.where(x < L/2, -x, ((x-1)**2)/2 - x)

>>> potential.allclose(analyticalArray, rtol = 2e-5, atol = 2e-5).value
1

and again view the result

>>> if __name__ == '__main__':
...     viewer.plot()
...     raw_input("Press any key to continue...")

Finally, we segregate all of the electrons to left side of the domain

.. math::

   C_{\text{e}^{-}} =
   \begin{cases}
       1& \text{for $x \le L/2$,} \\
       0& \text{for $x > L/2$.}
   \end{cases}

>>> interstitials[0].setValue(1.)
>>> interstitials[0].setValue(0., where=x > L / 2.)

and again solve for the electrostatic potential

>>> potential.equation.solve(var = potential)

which has the analytical solution

.. math::

   \psi(x) =
   \begin{cases}
       \frac{x^2}{2} - x& \text{for $x \le L/2$,} \\
       -\frac{1}{2}& \text{for $x > L/2$.}
   \end{cases}

We again verify that the correct equilibrium is attained

>>> analyticalArray = numerix.where(x < 1, (x**2)/2 - x, -0.5)

>>> potential.allclose(analyticalArray, rtol = 2e-5, atol = 2e-5).value
1

and again view the result

>>> if __name__ == '__main__':
...     viewer.plot()
"""
__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())

    raw_input("finished")
