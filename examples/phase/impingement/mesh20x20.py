#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 10/26/04 {9:00:00 PM} 
 #                                last update: 7/5/07 {8:21:37 PM}
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
 #  Description: 
 # 
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2004-10-26 JEG 1.0 original
 # ###################################################################
 ##

r"""

In the following examples, we solve the same set of equations as in::
    
    $ examples/phase/impingement/mesh40x1/input.py
    
with different initial conditions and a 2D mesh:

.. raw:: latex

   \IndexModule{parser}

..

    >>> from fipy.tools.parser import parse

    >>> numberOfElements = parse('--numberOfElements', action = 'store',
    ...                          type = 'int', default = 400)
    >>> numberOfSteps = parse('--numberOfSteps', action = 'store',
    ...                       type = 'int', default = 10)

.. raw:: latex

   \IndexFunction{sqrt}
   \IndexClass{Grid2D}

..
    
    >>> from fipy import *

    >>> steps = numberOfSteps
    >>> N = int(sqrt(numberOfElements))
    >>> L = 2.5 * N / 100.
    >>> dL = L / N
    >>> mesh = Grid2D(dx=dL, dy=dL, nx=N, ny=N)

The initial conditions are given by

.. raw:: latex

    $ \phi = 1 $ and
    
    $$ \theta = \begin{cases}
    \frac{2 \pi}{3} & \text{for $x^2 - y^2 < L / 2$,} \\
    \frac{-2 \pi}{3} & \text{for $(x-L)^2 - y^2 < L / 2$,} \\
    \frac{-2 \pi}{3}+0.3 & \text{for $x^2 - (y-L)^2 < L / 2$,} \\
    \frac{2 \pi}{3} & \text{for $(x-L)^2 - (y-L)^2 < L / 2$.}
    \end{cases} $$

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

.. raw:: latex

   \IndexClass{CellVariable}

..

    >>> phase = CellVariable(name='phase field', mesh=mesh)

The orientation is initialized to a uniform value to denote the
randomly oriented liquid phase

.. raw:: latex

   \IndexClass{ModularVariable}
   \IndexConstant{\pi}{pi}

..

    >>> theta = ModularVariable(
    ...     name='theta',
    ...     mesh=mesh,
    ...     value=-pi + 0.0001,
    ...     hasOld=1
    ...     )

Four different solid circular domains are created at each corner of
the domain with appropriate orientations

    >>> x, y = mesh.getCellCenters()
    >>> for a, b, thetaValue in ((0., 0.,  2. * pi / 3.), 
    ...                          (L, 0., -2. * pi / 3.), 
    ...                          (0., L, -2. * pi / 3. + 0.3), 
    ...                          (L, L,  2. * pi / 3.)):
    ...     segment = (x - a)**2 + (y - b)**2 < (L / 2.)**2
    ...     phase.setValue(1., where=segment)
    ...     theta.setValue(thetaValue, where=segment)

The `phase` equation is built in the following way. The source term is
linearized in the manner demonstrated in `examples.phase.simple.input`
(Kobayashi, semi-implicit). Here we use a function to build the equation,
so that it can be reused later.

.. raw:: latex

   \IndexClass{TransientTerm}
   \IndexClass{ExplicitDiffusionTerm}
   \IndexClass{ImplicitSourceTerm}

..

    >>> def buildPhaseEquation(phase, theta):
    ...
    ...     mPhiVar = phase - 0.5 + temperature * phase * (1 - phase)
    ...     thetaMag = theta.getOld().getGrad().getMag()
    ...     implicitSource = mPhiVar * (phase - (mPhiVar < 0))
    ...     implicitSource += (2 * s + epsilon**2 * thetaMag) * thetaMag
    ...
    ...     return TransientTerm(phaseTransientCoeff) == \
    ...               ExplicitDiffusionTerm(alpha**2) \
    ...               - ImplicitSourceTerm(implicitSource) \
    ...               + (mPhiVar > 0) * mPhiVar * phase

    >>> phaseEq = buildPhaseEquation(phase, theta)

The `theta` equation is built in the following way. The details for
this equation are fairly involved, see J.A. Warren *et al.*. The main
detail is that a source must be added to correct for the
discretization of `theta` on the circle.  The source term requires the
evaluation of the face gradient without the modular operators.

.. raw:: latex

   \IndexFunction{exp}

..

    >>> def buildThetaEquation(phase, theta):
    ...
    ...     phaseMod = phase + ( phase < thetaSmallValue ) * thetaSmallValue
    ...     phaseModSq = phaseMod * phaseMod
    ...     expo = epsilon * beta * theta.getGrad().getMag()
    ...     expo = (expo < 100.) * (expo - 100.) + 100.
    ...     pFunc = 1. + exp(-expo) * (mu / epsilon - 1.)
    ...
    ...     phaseFace = phase.getArithmeticFaceValue()
    ...     phaseSq = phaseFace * phaseFace
    ...     gradMag = theta.getFaceGrad().getMag()
    ...     eps = 1. / gamma / 10.
    ...     gradMag += (gradMag < eps) * eps
    ...     IGamma = (gradMag > 1. / gamma) * (1 / gradMag - gamma) + gamma
    ...     diffusionCoeff = phaseSq * (s * IGamma + epsilon**2)
    ...
    ...     thetaGradDiff = theta.getFaceGrad() - theta.getFaceGradNoMod()
    ...     sourceCoeff = (diffusionCoeff * thetaGradDiff).getDivergence()
    ...
    ...     return TransientTerm(thetaTransientCoeff * phaseModSq * pFunc) == \
    ...                ImplicitDiffusionTerm(diffusionCoeff) \
    ...                + sourceCoeff

.. raw:: latex

   \IndexClass{ImplicitDiffusionTerm}

..
    
    >>> thetaEq = buildThetaEquation(phase, theta)

If the example is run interactively, we create viewers for the phase
and orientation variables. Rather than viewing the raw orientation,
which is not meaningful in the liquid phase, we weight the orientation
by the phase

.. raw:: latex

   \IndexModule{viewers}

..

    >>> if __name__ == '__main__':
    ...     phaseViewer = viewers.make(vars=phase,
    ...                                limits={'datamin': 0., 'datamax': 1.})
    ...     thetaProd = -pi + phase * (theta + pi)
    ...     thetaProductViewer = make(vars=thetaProd,
    ...                               limits={'datamin': -pi, 
    ...                                       'datamax': pi})
    ...     phaseViewer.plot()
    ...     thetaProductViewer.plot()

The solution will be tested against data that was created with ``steps
= 10`` with a FORTRAN code written by Ryo Kobayashi for phase field
modeling. The following code opens the file `mesh20x20.gz` extracts the
data and compares it with the `theta` variable.

    >>> import os
    >>> import gzip
    >>> filestream = gzip.open(os.path.splitext(__file__)[0] + '.gz','r')
    >>> import cPickle
    >>> testData = cPickle.load(filestream)
    >>> filestream.close()
    
.. raw:: latex

   \IndexFunction{resize}

..

    >>> testData = resize(testData, (mesh.getNumberOfCells(),))
    
We step the solution in time, plotting as we go if running interactively,

    >>> for i in range(steps):
    ...     theta.updateOld()
    ...     phase.updateOld()
    ...     thetaEq.solve(theta, dt=timeStepDuration)
    ...     phaseEq.solve(phase, dt=timeStepDuration)
    ...     if __name__ == '__main__':
    ...         phaseViewer.plot()
    ...         thetaProductViewer.plot()
    
The solution is compared against Ryo Kobayashi's test data

    >>> print theta.allclose(testData, rtol=1e-7, atol=1e-7)
    1

The following code shows how to restart a simulation from some saved
data. First, reset the variables to their original values.

    >>> phase.setValue(0)
    >>> theta.setValue(-pi + 0.0001)
    >>> x, y = mesh.getCellCenters()
    >>> for a, b, thetaValue in ((0., 0.,  2. * pi / 3.), 
    ...                          (L, 0., -2. * pi / 3.), 
    ...                          (0., L, -2. * pi / 3. + 0.3), 
    ...                          (L, L,  2. * pi / 3.)):
    ...     segment = (x - a)**2 + (y - b)**2 < (L / 2.)**2
    ...     phase.setValue(1., where=segment)
    ...     theta.setValue(thetaValue, where=segment)

Step through half the time steps.

    >>> for i in range(steps / 2):
    ...     theta.updateOld()
    ...     phase.updateOld()
    ...     thetaEq.solve(theta, dt=timeStepDuration)
    ...     phaseEq.solve(phase, dt=timeStepDuration)

We confirm that the solution has not yet converged to that given by 
Ryo Kobayashi's FORTRAN code:

    >>> print theta.allclose(testData)
    0

We save the variables to disk.

.. raw:: latex

   \IndexModule{dump}

..

    >>> (f, filename) = dump.write({'phase' : phase, 'theta' : theta}, extension = '.gz')
    
and then recall them to test the data pickling mechanism

    >>> data = dump.read(filename, f)
    >>> newPhase = data['phase']
    >>> newTheta = data['theta']
    >>> newThetaEq = buildThetaEquation(newPhase, newTheta)
    >>> newPhaseEq = buildPhaseEquation(newPhase, newTheta)

and finish the iterations,

    >>> for i in range(steps / 2):
    ...     newTheta.updateOld()
    ...     newPhase.updateOld()
    ...     newThetaEq.solve(newTheta, dt=timeStepDuration)
    ...     newPhaseEq.solve(newPhase, dt=timeStepDuration)

The solution is compared against Ryo Kobayashi's test data

    >>> print newTheta.allclose(testData, rtol=1e-7)
    1
    
"""
__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())
    
    raw_input('finished')

