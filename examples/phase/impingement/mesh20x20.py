r"""Solve for the impingement of four grains in two dimensions.

In the following examples, we solve the same set of equations as in
:mod:`examples.phase.impingement.mesh40x1`
with different initial conditions and a 2D mesh:

.. index::
   pair: module; fipy.tools.parser

>>> from fipy.tools.parser import parse

>>> numberOfElements = parse('--numberOfElements', action = 'store',
...                          type = 'int', default = 400)
>>> numberOfSteps = parse('--numberOfSteps', action = 'store',
...                       type = 'int', default = 10)

.. index::
   single: sqrt
   single: Grid2D

>>> from fipy import CellVariable, ModularVariable, Grid2D, TransientTerm, DiffusionTerm, ExplicitDiffusionTerm, ImplicitSourceTerm, GeneralSolver, Viewer
>>> from fipy.tools import numerix, dump

>>> steps = numberOfSteps
>>> N = int(numerix.sqrt(numberOfElements))
>>> L = 2.5 * N / 100.
>>> dL = L / N
>>> mesh = Grid2D(dx=dL, dy=dL, nx=N, ny=N)

The initial conditions are given by
:math:`\phi = 1` and

.. math::

   \theta = \begin{cases}
   \frac{2 \pi}{3} & \text{for $x^2 - y^2 < L / 2$,} \\
   \frac{-2 \pi}{3} & \text{for $(x-L)^2 - y^2 < L / 2$,} \\
   \frac{-2 \pi}{3}+0.3 & \text{for $x^2 - (y-L)^2 < L / 2$,} \\
   \frac{2 \pi}{3} & \text{for $(x-L)^2 - (y-L)^2 < L / 2$.}
   \end{cases}

This defines four solid regions with different
orientations. Solidification occurs and then boundary wetting occurs
where the orientation varies.

The parameters for this example are

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

>>> temperature = 10.

and is initialized to liquid everywhere

.. index::
   single: CellVariable

>>> phase = CellVariable(name='phase field', mesh=mesh)

The orientation is initialized to a uniform value to denote the
randomly oriented liquid phase

.. index::
    single: ModularVariable
    single: :math:`\pi`
    single: pi

>>> theta = ModularVariable(
...     name='theta',
...     mesh=mesh,
...     value=-numerix.pi + 0.0001,
...     hasOld=1
...     )

Four different solid circular domains are created at each corner of
the domain with appropriate orientations

>>> x, y = mesh.cellCenters
>>> for a, b, thetaValue in ((0., 0.,  2. * numerix.pi / 3.),
...                          (L, 0., -2. * numerix.pi / 3.),
...                          (0., L, -2. * numerix.pi / 3. + 0.3),
...                          (L, L,  2. * numerix.pi / 3.)):
...     segment = (x - a)**2 + (y - b)**2 < (L / 2.)**2
...     phase.setValue(1., where=segment)
...     theta.setValue(thetaValue, where=segment)

The ``phase`` equation is built in the following way. The source term is
linearized in the manner demonstrated in :mod:`examples.phase.simple`
(Kobayashi, semi-implicit). Here we use a function to build the equation,
so that it can be reused later.

.. index::
   single: TransientTerm
   single: ExplicitDiffusionTerm
   single: ImplicitSourceTerm

>>> def buildPhaseEquation(phase, theta):
... 
...     mPhiVar = phase - 0.5 + temperature * phase * (1 - phase)
...     thetaMag = theta.old.grad.mag
...     implicitSource = mPhiVar * (phase - (mPhiVar < 0))
...     implicitSource += (2 * s + epsilon**2 * thetaMag) * thetaMag
... 
...     return TransientTerm(phaseTransientCoeff) == \
...               ExplicitDiffusionTerm(alpha**2) \
...               - ImplicitSourceTerm(implicitSource) \
...               + (mPhiVar > 0) * mPhiVar * phase

>>> phaseEq = buildPhaseEquation(phase, theta)

The ``theta`` equation is built in the following way. The details for
this equation are fairly involved, see J. A. Warren *et al.*. The main
detail is that a source must be added to correct for the
discretization of ``theta`` on the circle.  The source term requires the
evaluation of the face gradient without the modular operators.

.. index::
   single: exp

>>> def buildThetaEquation(phase, theta):
... 
...     phaseMod = phase + ( phase < thetaSmallValue ) * thetaSmallValue
...     phaseModSq = phaseMod * phaseMod
...     expo = epsilon * beta * theta.grad.mag
...     expo = (expo < 100.) * (expo - 100.) + 100.
...     pFunc = 1. + numerix.exp(-expo) * (mu / epsilon - 1.)
... 
...     phaseFace = phase.arithmeticFaceValue
...     phaseSq = phaseFace * phaseFace
...     gradMag = theta.faceGrad.mag
...     eps = 1. / gamma / 10.
...     gradMag += (gradMag < eps) * eps
...     IGamma = (gradMag > 1. / gamma) * (1 / gradMag - gamma) + gamma
...     diffusionCoeff = phaseSq * (s * IGamma + epsilon**2)
... 
...     thetaGradDiff = theta.faceGrad - theta.faceGradNoMod
...     sourceCoeff = (diffusionCoeff * thetaGradDiff).divergence
... 
...     return TransientTerm(thetaTransientCoeff * phaseModSq * pFunc) == \
...                DiffusionTerm(diffusionCoeff) \
...                + sourceCoeff

>>> thetaEq = buildThetaEquation(phase, theta)

If the example is run interactively, we create viewers for the phase
and orientation variables. Rather than viewing the raw orientation,
which is not meaningful in the liquid phase, we weight the orientation
by the phase

.. index::
   pair: module; fipy.viewers

>>> if __name__ == '__main__':
...     phaseViewer = Viewer(vars=phase, datamin=0., datamax=1.)
...     thetaProd = -numerix.pi + phase * (theta + numerix.pi)
...     thetaProductViewer = Viewer(vars=thetaProd,
...                                 datamin=-numerix.pi, datamax=numerix.pi)
...     phaseViewer.plot()
...     thetaProductViewer.plot()

The solution will be tested against data that was created with ``steps
= 10`` with a FORTRAN code written by Ryo Kobayashi for phase field
modeling. The following code opens the file :file:`mesh20x20.gz` extracts the
data and compares it with the `theta` variable.

.. index::
   single: loadtxt

>>> import os
>>> from future.utils import text_to_native_str
>>> testData = numerix.loadtxt(os.path.splitext(__file__)[0] + text_to_native_str('.gz')).flat

We step the solution in time, plotting as we go if running interactively,

>>> from builtins import range
>>> for i in range(steps):
...     theta.updateOld()
...     thetaEq.solve(theta, dt=timeStepDuration, solver=GeneralSolver(iterations=2000, tolerance=1e-15))
...     phaseEq.solve(phase, dt=timeStepDuration, solver=GeneralSolver(iterations=2000, tolerance=1e-15))
...     if __name__ == '__main__':
...         phaseViewer.plot()
...         thetaProductViewer.plot()

The solution is compared against Ryo Kobayashi's test data

>>> print(theta.allclose(testData, rtol=1e-7, atol=1e-7))
1

The following code shows how to restart a simulation from some saved
data. First, reset the variables to their original values.

>>> phase.setValue(0)
>>> theta.setValue(-numerix.pi + 0.0001)
>>> x, y = mesh.cellCenters
>>> for a, b, thetaValue in ((0., 0.,  2. * numerix.pi / 3.),
...                          (L, 0., -2. * numerix.pi / 3.),
...                          (0., L, -2. * numerix.pi / 3. + 0.3),
...                          (L, L,  2. * numerix.pi / 3.)):
...     segment = (x - a)**2 + (y - b)**2 < (L / 2.)**2
...     phase.setValue(1., where=segment)
...     theta.setValue(thetaValue, where=segment)

Step through half the time steps.

>>> from builtins import range
>>> for i in range(steps // 2):
...     theta.updateOld()
...     thetaEq.solve(theta, dt=timeStepDuration, solver=GeneralSolver(iterations=2000, tolerance=1e-15))
...     phaseEq.solve(phase, dt=timeStepDuration, solver=GeneralSolver(iterations=2000, tolerance=1e-15))

We confirm that the solution has not yet converged to that given by
Ryo Kobayashi's FORTRAN code:

>>> print(theta.allclose(testData))
0

We save the variables to disk.

.. index::
   pair: module; fipy.tools.dump

>>> (f, filename) = dump.write({'phase' : phase, 'theta' : theta}, extension = '.gz')

and then recall them to test the data pickling mechanism

>>> data = dump.read(filename, f)
>>> newPhase = data['phase']
>>> newTheta = data['theta']
>>> newThetaEq = buildThetaEquation(newPhase, newTheta)
>>> newPhaseEq = buildPhaseEquation(newPhase, newTheta)

and finish the iterations,

>>> from builtins import range
>>> for i in range(steps // 2):
...     newTheta.updateOld()
...     newThetaEq.solve(newTheta, dt=timeStepDuration, solver=GeneralSolver(iterations=2000, tolerance=1e-15))
...     newPhaseEq.solve(newPhase, dt=timeStepDuration, solver=GeneralSolver(iterations=2000, tolerance=1e-15))

The solution is compared against Ryo Kobayashi's test data

>>> print(newTheta.allclose(testData, rtol=1e-7))
1
"""
from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy import input

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())

    input('finished')
