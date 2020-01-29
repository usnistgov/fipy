r"""

In this example, a phase equation is solved in one dimension with a
misorientation between two solid domains.
The phase equation is given by:

.. math::

   \tau_{\phi} \frac{\partial \phi}{\partial t}
   = \alpha^2 \nabla^2 \phi + \phi ( 1 - \phi ) m_1 ( \phi , T)
   - 2 s \phi | \nabla \theta | - \epsilon^2 \phi | \nabla \theta |^2

where

.. math::

   m_1(\phi, T) = \phi - \frac{1}{2} - T \phi ( 1 - \phi )

The initial conditions are:

.. math::

   \phi &= 1 \qquad \forall x \\
   \theta &= \begin{cases}
   1 & \text{for $(x - L / 2)^2 + (y - L / 2)^2 > (L / 4)^2$} \\
   0 & \text{for $(x - L / 2)^2 + (y - L / 2)^2 \le (L / 4)^2$}
   \end{cases} \\
   T &= 1 \qquad \forall x

and boundary conditions
:math:`\phi = 1` for :math:`x = 0` and :math:`x = L`.

.. Further details of the numerical method for this problem can be found in
   "Extending Phase Field Models of Solidification to Polycrystalline
   Materials", J. A. Warren *et al.*, *Acta Materialia*, **51** (2003)
   6035-6058.

Here the phase equation is solved with an explicit technique.

The solution is allowed to evolve for ``steps = 100`` time steps.

>>> from builtins import range
>>> for step in range(steps):
...     phaseEq.solve(phase, dt = timeStepDuration)

The solution is compared with test data. The test data was created
with a FORTRAN code written by Ryo Kobayashi for phase field
modeling. The following code opens the file :file:`circle.gz` extracts the
data and compares it with the ``phase`` variable.

>>> import os
>>> from future.utils import text_to_native_str
>>> testData = numerix.loadtxt(os.path.splitext(__file__)[0] + text_to_native_str('.gz'))
>>> print(phase.allclose(testData))
1

"""
from __future__ import division
from __future__ import unicode_literals
from builtins import range
__docformat__ = 'restructuredtext'

from fipy import input
from fipy import CellVariable, ModularVariable, Grid2D, TransientTerm, ExplicitDiffusionTerm, ImplicitSourceTerm, Viewer
from fipy.tools import numerix

if __name__ == '__main__':
    steps = 100
else:
    steps = 10

timeStepDuration = 0.02
L = 1.5
nx = 100
ny = 100
temperature = 1.
phaseTransientCoeff = 0.1
epsilon = 0.008
s = 0.01
alpha = 0.015

dx = L / nx
dy = L / ny

mesh = Grid2D(dx, dy, nx, ny)

phase = CellVariable(name = 'PhaseField', mesh = mesh, value = 1.)

theta = ModularVariable(name = 'Theta', mesh = mesh, value = 1.)
x, y = mesh.cellCenters
theta.setValue(0., where=(x - L / 2.)**2 + (y - L / 2.)**2 < (L / 4.)**2)

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
    input('finished')
