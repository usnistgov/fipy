#!/usr/bin/env python

##
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "mesh1D.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # of Standards and Technology, an agency of the Federal Government.
 # Pursuant to title 17 section 105 of the United States Code,
 # United States Code this software is not subject to copyright
 # protection, and this software is considered to be in the public domain.
 # FiPy is an experimental system.
 # NIST assumes no responsibility whatsoever for its use by whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 #
 # To the extent that NIST may hold copyright in countries other than the
 # United States, you are hereby granted the non-exclusive irrevocable and
 # unconditional right to print, publish, prepare derivative works and
 # distribute this software, in any medium, or authorize others to do so on
 # your behalf, on a royalty-free basis throughout the world.
 #
 # You may improve, modify, and create derivative works of the software or
 # any portion of the software, and you may copy and distribute such
 # modifications or works.  Modified works should carry a notice stating
 # that you changed the software and should note the date and nature of any
 # such change.  Please explicitly acknowledge the National Institute of
 # Standards and Technology as the original source.
 #
 # This software can be redistributed and/or modified freely provided that
 # any derivative works bear some notice that they are derived from it, and
 # any modified versions bear some notice that they have been modified.
 # ========================================================================
 #
 # ###################################################################
 ##

r"""

In this example a phase equation is solved in 1 dimension with a
misorientation present. The phase equation is given by:

.. math::

   \tau_{\phi} \frac{\partial \phi}{\partial t}
   = \alpha^2 \nabla^2 \phi + \phi ( 1 - \phi ) m_1 ( \phi , T)
   - 2 s \phi | \nabla \theta | - \epsilon^2 \phi | \nabla \theta |^2

where

.. math::

   m_1(\phi, T) = \phi - \frac{1}{2} - T \phi ( 1 - \phi )

The initial conditions are:

.. math::

   \phi &= 1 \qquad \text{for $0 \le x \le L$} \\
   \theta &= \begin{cases}
   1 & \text{for $0 \le x \le L/2$} \\
   0 & \text{for $L/2 < x \le L$}
   \end{cases} \\
   T &= 1 \qquad \text{for $0 \le x \le L$}

and boundary conditions
:math:`\phi = 1` for :math:`x = 0` and :math:`x = L`.

.. Further details of the numerical method for this problem can be found in
   "Extending Phase Field Models of Solidification to Polycrystalline
   Materials", J.A. Warren *et al.*, *Acta Materialia*, **51** (2003)
   6035-6058.

Here the phase equation is solved with an explicit technique.

The solution is allowed to evolve for ``steps = 100`` time steps.

>>> for step in range(steps):
...     phaseEq.solve(phase, dt = timeStepDuration)

The solution is compared with test data. The test data was created
with a FORTRAN code written by Ryo Kobayashi for phase field
modeling. The following code opens the file :file:`mesh1D.gz` extracts the
data and compares it with the ``theta`` variable.

>>> import os
>>> testData = numerix.loadtxt(os.path.splitext(__file__)[0] + '.gz')
>>> print phase.allclose(testData)
1
"""
__docformat__ = 'restructuredtext'

from fipy import CellVariable, ModularVariable, Grid1D, TransientTerm, ExplicitDiffusionTerm, ImplicitSourceTerm, Viewer
from fipy.tools import numerix

steps = 100
timeStepDuration = 0.02
L = 1.5
nx = 100
temperature = 1.
phaseTransientCoeff = 0.1
epsilon = 0.008
s = 0.01
alpha = 0.015
temperature = 1.

dx = L / nx

mesh = Grid1D(dx = dx, nx = nx)

phase = CellVariable(name = 'PhaseField', mesh = mesh, value = 1.)

theta = ModularVariable(name = 'Theta', mesh = mesh, value = 1.)
theta.setValue(0., where=mesh.cellCenters[0] > L / 2.)

mPhiVar = phase - 0.5 + temperature * phase * (1 - phase)
thetaMag = theta.old.grad.mag
implicitSource = mPhiVar * (phase - (mPhiVar < 0))
implicitSource += (2 * s + epsilon**2 * thetaMag) * thetaMag

phaseEq = TransientTerm(phaseTransientCoeff) == \
          ExplicitDiffusionTerm(alpha**2) \
          - ImplicitSourceTerm(implicitSource) \
          + (mPhiVar > 0) * mPhiVar * phase

if __name__ == '__main__':

   phaseViewer = Viewer(vars = phase)
   phaseViewer.plot()
   for step in range(steps):
      phaseEq.solve(phase, dt = timeStepDuration)
      phaseViewer.plot()
   raw_input('finished')
