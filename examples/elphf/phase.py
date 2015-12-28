#!/usr/bin/env python

##
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "input.py"
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
A simple 1D example to test the setup of the phase field equation.

We rearrange Eq. :eq:`elphf:phase` to

.. math::

   \frac{1}{M_\xi}\frac{\partial \xi}{\partial t}
   &=
   \kappa_{\xi}\nabla^2 \xi
   +
   \frac{\epsilon'(\xi)}{2}\left(\nabla\phi\right)^2
   \\
   &\qquad -
   \left[
       p'(\xi) \Delta\mu_n^\circ
       + g'(\xi) W_n
   \right]
   -
   \sum_{j=2}^{n-1} C_j \left[
       p'(\xi) \Delta\mu_{jn}^\circ
       + g'(\xi) W_{jn}
   \right]
   -
   C_{\text{e}^{-}} \left[
       p'(\xi) \Delta\mu_{\text{e}^{-}}^\circ
       + g'(\xi) W_{\text{e}^{-}}
   \right]

The single-component phase field governing equation can be represented as

.. math::

   \frac{1}{M_\xi} \frac{\partial \xi}{\partial t}
   =  \kappa_\xi \nabla^2 \xi - 2\xi(1-\xi)(1-2\xi) W

where :math:`\xi` is the phase field, :math:`t` is time, :math:`M_\xi` is the
phase field mobility, :math:`\kappa_\xi` is the phase field gradient energy
coefficient, and :math:`W` is the phase field barrier energy.

We solve the problem on a 1D mesh

>>> from fipy import CellVariable, Grid1D, TransientTerm, DiffusionTerm, ImplicitSourceTerm, Viewer
>>> from fipy.tools import numerix

>>> nx = 400
>>> dx = 0.01
>>> L = nx * dx
>>> mesh = Grid1D(dx = dx, nx = nx)

We create the phase field

>>> phase = CellVariable(mesh = mesh, name = 'xi')
>>> phase.mobility = numerix.inf
>>> phase.gradientEnergy = 0.025

Although we are not interested in them for this problem, we create one field to
represent the "solvent" component (1 everywhere)

>>> class ComponentVariable(CellVariable):
...     def copy(self):
...         new = self.__class__(mesh = self.mesh,
...                              name = self.name,
...                              value = self.value)
...         new.standardPotential = self.standardPotential
...         new.barrier = self.barrier
...         return new

>>> solvent = ComponentVariable(mesh = mesh, name = 'Cn', value = 1.)
>>> solvent.standardPotential = 0.
>>> solvent.barrier = 1.

and one field to represent the electrostatic potential (0 everywhere)

>>> potential = CellVariable(mesh = mesh, name = 'phi', value = 0.)
>>> permittivityPrime = 0.

We'll have no substitutional species and no interstitial species in this first example

>>> substitutionals = []
>>> interstitials = []

>>> for component in substitutionals:
...     solvent -= component

>>> phase.equation = TransientTerm(coeff = 1/phase.mobility) \
...     == DiffusionTerm(coeff = phase.gradientEnergy) \
...     - (permittivityPrime / 2.) \
...        * potential.grad.dot(potential.grad)

>>> enthalpy = solvent.standardPotential
>>> barrier = solvent.barrier
>>> for component in substitutionals + interstitials:
...     enthalpy += component * component.standardPotential
...     barrier += component * component.barrier

We linearize the source term in the same way as in :mod:`example.phase.simple.input1D`.

>>> mXi = -(30 * phase * (1. - phase) * enthalpy \
...         +  4 * (0.5 - phase) * barrier)
>>> dmXidXi = (-60 * (0.5 - phase) * enthalpy + 4 * barrier)
>>> S1 = dmXidXi * phase * (1 - phase) + mXi * (1 - 2 * phase)
>>> S0 = mXi * phase * (1 - phase) - phase * S1

>>> phase.equation -= S0 + ImplicitSourceTerm(coeff = S1)

.. note::

   Adding a :class:`~fipy.terms.term.Term` to an equation formed with ``==`` will
   add to the left-hand side of the equation and subtracting a
   :class:`~fipy.terms.term.Term` will add to the right-hand side of the
   equation

We separate the phase field into electrode and electrolyte regimes

>>> phase.setValue(1.)
>>> phase.setValue(0., where=mesh.cellCenters[0] > L / 2)

Even though we are solving the steady-state problem (:math:`M_\phi = \infty`) we
still must sweep the solution several times to equilibrate

>>> for step in range(10):
...     phase.equation.solve(var = phase, dt=1.)

Since we have only a single component :math:`n`, with :math:`\Delta\mu_n^\circ =
0`, and the electrostatic potential is uniform, Eq. :eq:`elphf:phase` reduces to

.. math::

   \frac{1}{M_\xi}\frac{\partial \xi}{\partial t}
   = \kappa_{\xi}\nabla^2 \xi
   - g'(\xi) W_n

which we know from :mod:`examples.phase.simple.input1D` has the analytical
solution

.. math::

   \xi(x) = \frac{1}{2}(1 - \tanh\frac{x - L/2}{2d})

with an interfacial thickness :math:`d = \sqrt{\kappa_{\xi}/2W_n}`.

We verify that the correct equilibrium solution is attained

>>> x = mesh.cellCenters[0]

>>> d = numerix.sqrt(phase.gradientEnergy / (2 * solvent.barrier))
>>> analyticalArray = (1. - numerix.tanh((x - L/2.)/(2 * d))) / 2.

>>> phase.allclose(analyticalArray, rtol = 1e-4, atol = 1e-4).value
1


If we are running interactively, we plot the error

>>> if __name__ == '__main__':
...     viewer = Viewer(vars = (phase - \
...         CellVariable(name = "analytical", mesh = mesh,
...                      value = analyticalArray),))
...     viewer.plot()

.. image:: phase/error.*
   :scale: 50
   :align: center
   :alt: error in solution to steady-state phase field equation
"""
__docformat__ = 'restructuredtext'

if __name__ == '__main__':
##     from fipy.tools.profiler.profiler import Profiler
##     from fipy.tools.profiler.profiler import calibrate_profiler

    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())

##     fudge = calibrate_profiler(10000)
##     profile = Profiler('profile', fudge=fudge)

##     profile.stop()

    raw_input("finished")
