#!/usr/bin/env python

##
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "mesh40x1.py"
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

r"""Solve for the impingement of two grains in one dimension.

In this example we solve a coupled phase and orientation equation on a one
dimensional grid. This is another aspect of the model of Warren, Kobayashi,
Lobkovsky and Carter :cite:`WarrenPolycrystal`

.. index:: Grid1D

>>> from fipy import CellVariable, ModularVariable, Grid1D, TransientTerm, DiffusionTerm, ExplicitDiffusionTerm, ImplicitSourceTerm, GeneralSolver, Viewer
>>> from fipy.tools import numerix

>>> nx = 40
>>> Lx = 2.5 * nx / 100.
>>> dx = Lx / nx
>>> mesh = Grid1D(dx=dx, nx=nx)

This problem simulates the wet boundary that forms between grains of different
orientations. The phase equation is given by

.. math::

   \tau_{\phi} \frac{\partial \phi}{\partial t}
   = \alpha^2 \nabla^2 \phi + \phi ( 1 - \phi ) m_1 ( \phi , T)
   - 2 s \phi | \nabla \theta | - \epsilon^2 \phi | \nabla \theta |^2

where

.. math::

   m_1(\phi, T) = \phi - \frac{1}{2} - T \phi ( 1 - \phi )

and the orientation equation is given by

.. math::

   P(\epsilon | \nabla \theta |) \tau_{\theta} \phi^2
   \frac{\partial \theta}{\partial t}
   = \nabla \cdot \left[ \phi^2 \left( \frac{s}{| \nabla \theta |}
   + \epsilon^2 \right) \nabla \theta \right]

where

.. math::

   P(w) = 1 - \exp{(-\beta w)} + \frac{\mu}{\epsilon} \exp{(-\beta w)}

The initial conditions for this problem are set such that
:math:`\phi = 1` for :math:`0 \le x \le L_x` and

.. math::

   \theta = \begin{cases}
   1 & \text{for $0 \le x < L_x / 2$,} \\
   0 & \text{for $L_x / 2 \le x \le L_x$.}
   \end{cases}

.. Further details of the numerical method for this problem can be found in
   "Extending Phase Field Models of Solidification to Polycrystalline
   Materials", J.A. Warren *et al.*, *Acta Materialia*, **51** (2003)
   6035-6058.

Here the phase and orientation equations are solved with an
explicit and implicit technique respectively.

The parameters for these equations are

>>> timeStepDuration = 0.02
>>> phaseTransientCoeff = 0.1
>>> thetaSmallValue = 1e-6
>>> beta = 1e5
>>> mu = 1e3
>>> thetaTransientCoeff = 0.01
>>> gamma= 1e3
>>> epsilon = 0.008
>>> s = 0.01
>>> alpha = 0.015

The system is held isothermal at

    >>> temperature = 1.

and is initially solid everywhere

.. index:: CellVariable

>>> phase = CellVariable(
...     name='phase field',
...     mesh=mesh,
...     value=1.
...     )

Because ``theta``
is an :math:`S^1`-valued variable (i.e. it maps to the circle) and thus
intrinsically has :math:`2\pi`-peridocity,
we must use :class:`~fipy.variables.modularVariable.ModularVariable` instead of a :class:`~fipy.variables.cellVariable.CellVariable`. A
:class:`~fipy.variables.modularVariable.ModularVariable` confines ``theta`` to
:math:`-\pi < \theta \le \pi` by adding or subtracting :math:`2\pi` where
necessary and by defining a new
subtraction operator between two angles.

>>> theta = ModularVariable(
...     name='theta',
...     mesh=mesh,
...     value=1.,
...     hasOld=1
...     )

The left and right halves of the domain are given different orientations.

>>> theta.setValue(0., where=mesh.cellCenters[0] > Lx / 2.)

The ``phase`` equation is built in the following way.

.. index:: TransientTerm, ExplicitDiffusionTerm, ImplicitSourceTerm

>>> mPhiVar = phase - 0.5 + temperature * phase * (1 - phase)

The source term is linearized in the manner demonstrated in
:mod:`examples.phase.simple` (Kobayashi, semi-implicit).

>>> thetaMag = theta.old.grad.mag
>>> implicitSource = mPhiVar * (phase - (mPhiVar < 0))
>>> implicitSource += (2 * s + epsilon**2 * thetaMag) * thetaMag

The ``phase`` equation is constructed.

>>> phaseEq = TransientTerm(phaseTransientCoeff) \
...   == ExplicitDiffusionTerm(alpha**2) \
...      - ImplicitSourceTerm(implicitSource) \
...      + (mPhiVar > 0) * mPhiVar * phase

The ``theta`` equation is built in the following way. The details for
this equation are fairly involved, see J.A. Warren *et al.*. The main
detail is that a source must be added to correct for the
discretization of ``theta`` on the circle.

.. index:: exp

>>> phaseMod = phase + ( phase < thetaSmallValue ) * thetaSmallValue
>>> phaseModSq = phaseMod * phaseMod
>>> expo = epsilon * beta * theta.grad.mag
>>> expo = (expo < 100.) * (expo - 100.) + 100.
>>> pFunc = 1. + numerix.exp(-expo) * (mu / epsilon - 1.)

>>> phaseFace = phase.arithmeticFaceValue
>>> phaseSq = phaseFace * phaseFace
>>> gradMag = theta.faceGrad.mag
>>> eps = 1. / gamma / 10.
>>> gradMag += (gradMag < eps) * eps
>>> IGamma = (gradMag > 1. / gamma) * (1 / gradMag - gamma) + gamma
>>> diffusionCoeff = phaseSq * (s * IGamma + epsilon**2)

The source term requires the evaluation of the face gradient without
the modular operator. ``theta``:meth:`~fipy.variables.modularVariable.ModularVariable.getFaceGradNoMod`
evelautes the gradient without modular arithmetic.

>>> thetaGradDiff = theta.faceGrad - theta.faceGradNoMod
>>> sourceCoeff = (diffusionCoeff * thetaGradDiff).divergence

Finally the ``theta`` equation can be constructed.

>>> thetaEq = TransientTerm(thetaTransientCoeff * phaseModSq * pFunc) == \
...           DiffusionTerm(diffusionCoeff) \
...           + sourceCoeff

If the example is run interactively, we create viewers for the phase
and orientation variables.

.. index::
   module: fipy.viewers

.. index:: :math:`\pi`, pi

>>> if __name__ == '__main__':
...     phaseViewer = Viewer(vars=phase, datamin=0., datamax=1.)
...     thetaProductViewer = Viewer(vars=theta,
...                                 datamin=-numerix.pi, datamax=numerix.pi)
...     phaseViewer.plot()
...     thetaProductViewer.plot()

we iterate the solution in time, plotting as we go if running interactively,

>>> steps = 10
>>> for i in range(steps):
...     theta.updateOld()
...     thetaEq.solve(theta, dt = timeStepDuration)
...     phaseEq.solve(phase, dt = timeStepDuration)
...     if __name__ == '__main__':
...         phaseViewer.plot()
...         thetaProductViewer.plot()

The solution is compared with test data. The test data was created
with ``steps = 10`` with a FORTRAN code written by Ryo Kobayashi for
phase field modeling. The following code opens the file :file:`mesh40x1.gz`
extracts the data and compares it with the ``theta`` variable.

.. index:: loadtxt

>>> import os
>>> testData = numerix.loadtxt(os.path.splitext(__file__)[0] + '.gz')
>>> testData = CellVariable(mesh=mesh, value=testData)
>>> print theta.allclose(testData)
1
"""
__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())

    raw_input('finished')
