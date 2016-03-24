#!/usr/bin/env python

##
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "input2D.py"
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

r"""
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
   = \nabla\cdot D \nabla\left( \frac{\partial f}{\partial \phi}   - \epsilon^2 \nabla^2 \phi\right).

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

>>> from fipy import CellVariable, Grid2D, GaussianNoiseVariable, TransientTerm, DiffusionTerm, ImplicitSourceTerm, LinearLUSolver, Viewer
>>> from fipy.tools import numerix

(Note that all of the functionality of NumPy is imported along with :term:`FiPy`, although
much is augmented for :term:`FiPy`\'s needs.)

>>> if __name__ == "__main__":
...     nx = ny = 1000
... else:
...     nx = ny = 20
>>> mesh = Grid2D(nx=nx, ny=ny, dx=0.25, dy=0.25)
>>> phi = CellVariable(name=r"$\phi$", mesh=mesh)

We start the problem with random fluctuations about :math:`\phi = 1/2`

>>> phi.setValue(GaussianNoiseVariable(mesh=mesh,
...                                    mean=0.5,
...                                    variance=0.01))

:term:`FiPy` doesn't plot or output anything unless you tell it to:

>>> if __name__ == "__main__":
...     viewer = Viewer(vars=(phi,), datamin=0., datamax=1.)

For :term:`FiPy`, we need to perform the partial derivative
:math:`\partial f/\partial \phi`
manually and then put the equation in the canonical
form by decomposing the spatial derivatives
so that each :class:`~fipy.terms.term.Term` is of a single, even order:

.. math::

   \frac{\partial \phi}{\partial t}
    = \nabla\cdot D a^2 \left[ 1 - 6 \phi \left(1 - \phi\right)\right] \nabla \phi- \nabla\cdot D \nabla \epsilon^2 \nabla^2 \phi.

:term:`FiPy` would automatically interpolate
``D * a**2 * (1 - 6 * phi * (1 - phi))``
onto the faces, where the diffusive flux is calculated, but we obtain
somewhat more accurate results by performing a linear interpolation from
``phi`` at cell centers to ``PHI`` at face centers.
Some problems benefit from non-linear interpolations, such as harmonic or
geometric means, and :term:`FiPy` makes it easy to obtain these, too.

>>> PHI = phi.arithmeticFaceValue
>>> D = a = epsilon = 1.
>>> eq = (TransientTerm()
...       == DiffusionTerm(coeff=D * a**2 * (1 - 6 * PHI * (1 - PHI)))
...       - DiffusionTerm(coeff=(D, epsilon**2)))

Because the evolution of a spinodal microstructure slows with time, we
use exponentially increasing time steps to keep the simulation
"interesting". The :term:`FiPy` user always has direct control over the
evolution of their problem.

>>> dexp = -5
>>> elapsed = 0.
>>> if __name__ == "__main__":
...     duration = 1000.
... else:
...     duration = 1000.

>>> while elapsed < duration:
...     dt = min(100, numerix.exp(dexp))
...     elapsed += dt
...     dexp += 0.01
...     eq.solve(phi, dt=dt, solver=LinearLUSolver())
...     if __name__ == "__main__":
...         viewer.plot()
...     elif (max(phi.globalValue) > 0.7) and (min(phi.globalValue) < 0.3) and elapsed > 10.:
...         break

>>> print (max(phi.globalValue) > 0.7) and (min(phi.globalValue) < 0.3)
True

.. image:: mesh2D.*
   :width: 90%
   :align: center
   :alt: evolution of Cahn-Hilliard phase separation at t = 30, 100 and 1000

"""
__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())

    input('finished')
