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
.. attention::

   This example remains only for exact comparison against Ryo Kobayashi's
   FORTRAN code. See :mod:`examples.phase.anisotropy` for a better, although not
   numerically identical implementation.

In this example we solve a coupled phase and temperature equation to model
solidification, and eventually dendritic growth, based on the work of
Warren, Kobayashi, Lobkovsky and Carter :cite:`WarrenPolycrystal`.

We start from a circular seed in a 2D mesh:

.. index:: Grid2D

>>> from fipy import CellVariable, Grid2D, TransientTerm, DiffusionTerm, ExplicitDiffusionTerm, ImplicitSourceTerm, Viewer
>>> from fipy.tools import numerix

>>> numberOfCells = 40
>>> Length = numberOfCells * 2.5 / 100.
>>> nx = numberOfCells
>>> ny = numberOfCells
>>> dx = Length / nx
>>> dy = Length / ny
>>> radius = Length / 4.
>>> seedCenter = (Length / 2., Length / 2.)
>>> initialTemperature = -0.4
>>> mesh = Grid2D(dx=dx, dy=dy, nx=nx, ny=ny)

Dendritic growth will not be observed with this small test system. If
you wish to see dendritic growth reset the following parameters such
that ``numberOfCells = 500``, ``steps = 10000``, ``radius = dx * 5.``
``seedCenter = (0. , 0.)`` and ``initialTemperature = -0.5``.

The governing equation for the phase field is given by:

.. math::

   \tau_{\phi} \frac{\partial \phi}{\partial t}
   = \nabla \cdot \left[ D \nabla \phi + A \nabla \xi \right] +
   \phi ( 1 - \phi ) m ( \phi , T)

where

.. math::

   m(\phi, T)
   = \phi - \frac{1}{2} - \frac{ \kappa_1 }{ \pi } \arctan \left( \kappa_2 T \right).

The coefficients :math:`D` and :math:`A` are given by,

.. math::

   D = \alpha^2 \left[ 1 + c \beta \right]^2

and

.. math::

   A = \alpha^2 c \left[ 1 + c \beta \right] \beta_\psi

where :math:`\beta = \frac{ 1 - \Phi^2 } { 1 + \Phi^2}`,
:math:`\Phi = \tan \left( \frac{ N } { 2 } \psi \right) `,
:math:`\psi = \theta + \arctan \frac{ \phi_y } { \phi_x }` and
:math:`\xi_x = -\phi_y` and :math:`\xi_y = \phi_x`.

The governing equation for temperature is given by:

.. math::

   \frac{\partial T}{\partial t} = D_T \nabla^2 T + \frac{\partial \phi}{\partial t}

..  Further details of the numerical method for this problem can be found in
    "Extending Phase Field Models of Solidification to Polycrystalline
    Materials", J.A. Warren *et al.*, *Acta Materialia*, **51** (2003) 6035-6058.

Here the phase and temperature equations are solved with an explicit
and implicit technique, respectively.

The parameters for these equations are

>>> timeStepDuration = 5e-5
>>> tau = 3e-4
>>> alpha = 0.015
>>> c = 0.02
>>> N = 4.
>>> kappa1 = 0.9
>>> kappa2 = 20.
>>> tempDiffusionCoeff = 2.25
>>> theta = 0.

The ``phase`` variable is 0 for a liquid and ``1 for a solid.  Here,
the ``phase`` variable is initialized as a liquid,

.. index:: CellVariable

>>> phase = CellVariable(name='phase field', mesh=mesh, hasOld=1)

The ``hasOld`` flag keeps the old value of the variable. This is
necessary for a transient solution. In this example we wish to set up
an interior region that is solid.
The domain is seeded with a circular solidified region with parameters
``seedCenter`` and ``radius`` representing the center and radius of the
seed.

>>> x, y = mesh.cellCenters
>>> phase.setValue(1., where=((x - seedCenter[0])**2
...                           + (y - seedCenter[1])**2) < radius**2)

The temperature field is initialized to a value of -0.4 throughout:

>>> temperature = CellVariable(
...     name='temperature',
...     mesh=mesh,
...     value=initialTemperature,
...     hasOld=1)

The :math:`m(\phi, T)` variable

is created from the ``phase`` and ``temperature`` variables.

.. index:: :math:`\pi`, pi, arctan, arctan2, tan

>>> mVar = phase - 0.5 - kappa1 / numerix.pi * numerix.arctan(kappa2 * temperature)

The following section of code builds up the :math:`A` and :math:`D` coefficients.

>>> phaseY = phase.faceGrad.dot((0, 1))
>>> phaseX = phase.faceGrad.dot((1, 0))
>>> psi = theta + numerix.arctan2(phaseY, phaseX)
>>> Phi = numerix.tan(N * psi / 2)
>>> PhiSq = Phi**2
>>> beta = (1. - PhiSq) / (1. + PhiSq)
>>> betaPsi = -N * 2 * Phi / (1 + PhiSq)
>>> A = alpha**2 * c * (1.+ c * beta) * betaPsi
>>> D = alpha**2 * (1.+ c * beta)**2

The :math:`\nabla \xi` variable (``dxi``),
given by :math:`(\xi_x, \xi_y) = (-\phi_y, \phi_x)`,
is constructed by first obtaining :math:`\nabla \phi`

using :meth:`getFaceGrad`. The axes are rotated ninety degrees.

>>> dxi = phase.faceGrad.dot(((0, 1),(-1,0)))
>>> anisotropySource = (A * dxi).divergence

The phase equation can now be constructed.

.. index:: TransientTerm, ExplicitDiffusionTerm, ImplicitSourceTerm

>>> phaseEq = TransientTerm(tau) == ExplicitDiffusionTerm(D) + \
...     ImplicitSourceTerm(mVar * ((mVar < 0) - phase)) + \
...     ((mVar > 0.) * mVar * phase + anisotropySource)

The temperature equation is built in the following way,

>>> temperatureEq = TransientTerm() == \
...                 DiffusionTerm(tempDiffusionCoeff) + \
...                 (phase - phase.old) / timeStepDuration

If we are running this example interactively, we create viewers for
the phase and temperature fields

.. index::
   module: fipy.viewers

>>> if __name__ == '__main__':
...     phaseViewer = Viewer(vars=phase)
...     temperatureViewer = Viewer(vars=temperature,
...                                datamin=-0.5, datamax=0.5)
...     phaseViewer.plot()
...     temperatureViewer.plot()

we iterate the solution in time, plotting as we go if running interactively,

>>> steps = 10
>>> for i in range(steps):
...     phase.updateOld()
...     temperature.updateOld()
...     phaseEq.solve(phase, dt=timeStepDuration)
...     temperatureEq.solve(temperature, dt=timeStepDuration)
...     if i%10 == 0 and __name__ == '__main__':
...         phaseViewer.plot()
...         temperatureViewer.plot()

The solution is compared with test data. The test data was created for
``steps = 10`` with a FORTRAN code written by Ryo Kobayashi for phase
field modeling. The following code opens the file :file:`anisotropy.gz` extracts
the data and compares it with the `phase` variable.

.. index:: loadtxt, allclose

>>> import os
>>> testData = numerix.loadtxt(os.path.splitext(__file__)[0] + '.gz')
>>> print phase.allclose(testData)
1
"""
__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())

    raw_input('finished')
