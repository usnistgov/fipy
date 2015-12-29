#!/usr/bin/env python

##
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "mesh2DCoupled.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #    mail: NIST
 #     www: http://ctcms.nist.gov
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

r"""Solve the Cahn-Hilliard problem in two dimensions.

The spinodal decomposition phenomenon is a spontaneous separation of
an initially homogenous mixture into two distinct regions of different
properties (spin-up/spin-down, component A/component B). It is a
"barrierless" phase separation process, such that under the right
thermodynamic conditions, any fluctuation, no matter how small, will
tend to grow. This is in contrast to nucleation, where a fluctuation
must exceed some critical magnitude before it will survive and grow.
Spinodal decomposition can be described by the "Cahn-Hilliard"
equation (also known as
"conserved Ginsberg-Landau" or "model B" of Hohenberg & Halperin)

.. math::

   \frac{\partial \phi}{\partial t}
   = \nabla \cdot D \nabla \left( \frac{\partial f}{\partial \phi}   - \epsilon^2 \nabla^2 \phi \right).

where :math:`\phi` is a conserved order parameter, possibly representing
alloy composition or spin.
The double-well free energy function :math:`f = (a^2/2) \phi^2 (1 - \phi)^2`
penalizes states with intermediate values of :math:`\phi`
between 0 and 1. The gradient energy term :math:`\epsilon^2 \nabla^2\phi`,
on the other hand, penalizes sharp changes of :math:`\phi`.
These two competing effects result in the segregation
of :math:`\phi` into domains of 0 and 1, separated by abrupt, but
smooth, transitions. The parameters :math:`a` and :math:`\epsilon` determine the relative
weighting of the two effects and :math:`D` is a rate constant.

We can simulate this process in :term:`FiPy` with a simple script:

>>> from fipy import CellVariable, Grid2D, GaussianNoiseVariable, DiffusionTerm, TransientTerm, ImplicitSourceTerm, Viewer
>>> from fipy.tools import numerix

(Note that all of the functionality of NumPy is imported along with :term:`FiPy`, although
much is augmented for :term:`FiPy`\'s needs.)

>>> if __name__ == "__main__":
...     nx = ny = 20
... else:
...     nx = ny = 10
>>> mesh = Grid2D(nx=nx, ny=ny, dx=0.25, dy=0.25)
>>> phi = CellVariable(name=r"$\phi$", mesh=mesh)
>>> psi = CellVariable(name=r"$\psi$", mesh=mesh)

We start the problem with random fluctuations about :math:`\phi = 1/2`

>>> noise = GaussianNoiseVariable(mesh=mesh,
...                               mean=0.5,
...                               variance=0.01).value

>>> phi[:] = noise

:term:`FiPy` doesn't plot or output anything unless you tell it to:

>>> if __name__ == "__main__":
...     viewer = Viewer(vars=(phi,), datamin=0., datamax=1.)

We factor the Cahn-Hilliard equation into two 2nd-order PDEs and place
them in canonical form for :term:`FiPy` to solve them as a coupled set
of equations.

.. math::

   \frac{\partial \phi}{\partial t} &= \nabla\cdot D \nabla \psi \\
   \psi &= \left(\frac{\partial f}{\partial \phi} - \frac{\partial^2 f}{\partial \phi^2}\phi\right)_{\text{old}}
           + \frac{\partial^2 f}{\partial \phi^2}\phi - \epsilon^2 \nabla^2 \phi

The source term in :math:`\psi`, :math:`\frac{\partial f}{\partial \phi}`, is
expressed in linearized form after Taylor expansion at :math:`\phi=\phi_{\text{old}}`,
for the same reasons discussed in :mod:`examples.phase.simple`.
We need to perform the partial derivatives

.. math::

   \frac{\partial f}{\partial \phi} &= (a^2/2) 2 \phi (1 - \phi) (1 - 2 \phi) \\
   \frac{\partial^2 f}{\partial \phi^2} &= (a^2/2) 2 \left[1 - 6 \phi(1 - \phi)\right]

manually.

>>> D = a = epsilon = 1.
>>> dfdphi = a**2 * phi * (1 - phi) * (1 - 2 * phi)
>>> dfdphi_ = a**2 * (1 - phi) * (1 - 2 * phi)
>>> d2fdphi2 = a**2 * (1 - 6 * phi * (1 - phi))
>>> eq1 = (TransientTerm(var=phi) == DiffusionTerm(coeff=D, var=psi))
>>> eq2 = (ImplicitSourceTerm(coeff=1., var=psi)
...        == ImplicitSourceTerm(coeff=d2fdphi2, var=phi) - d2fdphi2 * phi + dfdphi
...        - DiffusionTerm(coeff=epsilon**2, var=phi))
>>> eq3 = (ImplicitSourceTerm(coeff=1., var=psi)
...        == ImplicitSourceTerm(coeff=dfdphi_, var=phi)
...        - DiffusionTerm(coeff=epsilon**2, var=phi))

>>> eq = eq1 & eq2

Because the evolution of a spinodal microstructure slows with time, we
use exponentially increasing time steps to keep the simulation
"interesting". The :term:`FiPy` user always has direct control over the
evolution of their problem.

>>> dexp = -5
>>> elapsed = 0.
>>> if __name__ == "__main__":
...     duration = 100.
... else:
...     duration = .5e-1

>>> while elapsed < duration:
...     dt = min(100, numerix.exp(dexp))
...     elapsed += dt
...     dexp += 0.01
...     eq.solve(dt=dt)
...     if __name__ == "__main__":
...         viewer.plot()

>>> if __name__ == '__main__':
...     raw_input("Coupled equations. Press <return> to proceed...")

.. image:: mesh2DCoupled.*
   :width: 90%
   :align: center

-----

These equations can also be solved in :term:`FiPy` using a vector
equation. The variables :math:`\phi` and :math:`\psi` are now stored in
a single variable

>>> var = CellVariable(mesh=mesh, elementshape=(2,))
>>> var[0] = noise

>>> if __name__ == "__main__":
...     viewer = Viewer(name=r"$\phi$", vars=var[0,], datamin=0., datamax=1.)

>>> D = a = epsilon = 1.
>>> v0 = var[0]
>>> dfdphi = a**2 * v0 * (1 - v0) * (1 - 2 * v0)
>>> dfdphi_ = a**2 * (1 - v0) * (1 - 2 * v0)
>>> d2fdphi2 = a**2 * (1 - 6 * v0 * (1 - v0))

The source terms have to be shaped correctly for a vector. The implicit source
coefficient has to have a shape of `(2, 2)` while the explicit source
has a shape `(2,)`

>>> source = (- d2fdphi2 * v0 + dfdphi) * (0, 1)
>>> impCoeff = d2fdphi2 * ((0, 0),
...                        (1., 0)) + ((0, 0),
...                                    (0, -1.))

This is the same equation as the previous definition of `eq`, but now in
a vector format.

>>> eq = (TransientTerm(((1., 0.),
...                     (0., 0.))) == DiffusionTerm([((0.,          D),
...                                                   (-epsilon**2, 0.))])
...                                   + ImplicitSourceTerm(impCoeff) + source)

>>> dexp = -5
>>> elapsed = 0.

>>> while elapsed < duration:
...     dt = min(100, numerix.exp(dexp))
...     elapsed += dt
...     dexp += 0.01
...     eq.solve(var=var, dt=dt)
...     if __name__ == "__main__":
...         viewer.plot()

>>> print numerix.allclose(var, (phi, psi))
True

"""
__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())

    raw_input('finished')
