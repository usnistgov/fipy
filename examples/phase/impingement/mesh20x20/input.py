#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 10/26/04 {9:00:00 PM} 
 #                                last update: 3/10/05 {5:12:29 PM}
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

    >>> from fipy.tools.parser import parse

    >>> numberOfElements = parse('--numberOfElements', action = 'store',
    ...                          type = 'int', default = 400)
    >>> numberOfSteps = parse('--numberOfSteps', action = 'store',
    ...                       type = 'int', default = 10)

    >>> steps = numberOfSteps
    >>> import Numeric
    >>> nx = int(Numeric.sqrt(numberOfElements))
    >>> ny = nx
    >>> Lx = 2.5 * nx / 100.
    >>> dx = Lx / nx
    >>> from fipy.meshes.grid2D import Grid2D
    >>> mesh = Grid2D(dx,dx,nx,nx)

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

    >>> from fipy.variables.cellVariable import CellVariable
    >>> phase = CellVariable(
    ...     name = 'PhaseField',
    ...     mesh = mesh,
    ...     value = 0.
    ...     )

The orientation is initialized to a uniform value to denote the
randomly oriented liquid phase

    >>> from Numeric import pi
    >>> from fipy.models.phase.theta.modularVariable import ModularVariable
    >>> theta = ModularVariable(
    ...     name = 'Theta',
    ...     mesh = mesh,
    ...     value = -pi + 0.0001,
    ...     hasOld = 1
    ...     )

Four different solid circular domains are created at each corner of
the domain with appropriate orientations

    >>> def cornerCircle(cell, a = 1., b = 1., Lx = Lx):
    ...     x = cell.getCenter()[0]
    ...     y = cell.getCenter()[1]
    ...     if ((x - a)**2 + (y - b)**2) < (Lx / 2.)**2:
    ...         return 1
    ...     else:
    ...         return 0
    
    >>> for a, b, thetaValue in ((0., 0.,  2. * pi / 3.), 
    ...                          (Lx, 0., -2. * pi / 3.), 
    ...                          (0., Lx, -2. * pi / 3. + 0.3), 
    ...                          (Lx, Lx,  2. * pi / 3.)):
    ...     cells = mesh.getCells(filter = cornerCircle, a = a, b = b)
    ...     phase.setValue(1., cells)
    ...     theta.setValue(thetaValue, cells)

The `phase` equation is built in the following way.

    >>> from fipy.terms.transientTerm import TransientTerm
    >>> from fipy.terms.explicitDiffusionTerm import ExplicitDiffusionTerm
    >>> from fipy.terms.implicitSourceTerm import ImplicitSourceTerm

    >>> mPhiVar = phase - 0.5 + temperature * phase * (1 - phase)

The source term is linearized in the manner demonstrated in
`examples.phase.simple.input` (Kobayashi, semi-implicit).

    >>> thetaMag = theta.getOld().getGrad().getMag()
    >>> implicitSource = mPhiVar * (phase - (mPhiVar < 0))
    >>> implicitSource += (2 * s + epsilon**2 * thetaMag) * thetaMag
    >>> phaseEq = TransientTerm(phaseTransientCoeff) - ExplicitDiffusionTerm(alpha**2)
    >>> phaseEq += ImplicitSourceTerm(implicitSource) - (mPhiVar > 0) * mPhiVar * phase

The `theta` equation is built in the following way. The details for
this equation are fairly involved, see J.A. Warren *et al.*. The main
detail is that a source must be added to correct for the
discretization of `theta` on the circle.

    >>> phaseMod = phase + ( phase < thetaSmallValue ) * thetaSmallValue
    >>> phaseModSq = phaseMod * phaseMod
    >>> expo = epsilon * beta * theta.getGrad().getMag()
    >>> expo = (expo < 100.) * (expo - 100.) + 100.
    >>> import fipy.tools.array
    >>> pFunc = 1. + fipy.tools.array.exp(-expo) * (mu / epsilon - 1.)

    >>> phaseFace = phase.getArithmeticFaceValue()
    >>> phaseSq = phaseFace * phaseFace
    >>> gradMag = theta.getFaceGrad().getMag()
    >>> eps = 1. / gamma / 10.
    >>> gradMag += (gradMag < eps) * eps
    >>> IGamma = (gradMag > 1. / gamma) * (1 / gradMag - gamma) + gamma
    >>> diffusionCoeff = phaseSq * (s * IGamma + epsilon**2)

The source term requires the evaluation of the face gradient without
the modular operators. Thus a new subclass of `CellVariable` is
created that uses the value of the `ModularVariable` but does not use
its operators.
    
    >>> class NonModularTheta(CellVariable):
    ...     def __init__(self, modVar):
    ...         CellVariable.__init__(self, mesh = modVar.getMesh())
    ...         self.modVar = self._requires(modVar)
    ...     def _calcValue(self):
    ...         self.value = self.modVar[:]
        
    >>> thetaNoMod = NonModularTheta(theta)
    >>> thetaGradDiff = theta.getFaceGrad() - thetaNoMod.getFaceGrad()
    >>> from fipy.models.phase.phase.addOverFacesVariable import AddOverFacesVariable
    >>> sourceCoeff = AddOverFacesVariable(faceGradient = thetaGradDiff, faceVariable = diffusionCoeff)

    >>> transientTerm = TransientTerm(thetaTransientCoeff * phaseModSq * pFunc)
    >>> from fipy.terms.implicitDiffusionTerm import ImplicitDiffusionTerm
    >>> diffusionTerm = ImplicitDiffusionTerm(diffusionCoeff)

    >>> thetaEq = transientTerm - diffusionTerm - sourceCoeff

If the example is run interactively, we create viewers for the phase
and orientation variables. Rather than viewing the raw orientation,
which is not meaningful in the liquid phase, we weight the orientation
by the phase

    >>> if __name__ == '__main__':
    ...     import fipy.viewers
    ...     phaseViewer = fipy.viewers.make(vars = phase, 
    ...                                     limits = {'datamin': 0., 'datamax': 1.})
    ...     from Numeric import pi
    ...     thetaProd = -pi + phase * (theta + pi)
    ...     thetaProductViewer = fipy.viewers.make(vars = thetaProd ,
    ...                                            limits = {'datamin': -pi, 'datamax': pi})
    ...     phaseViewer.plot()
    ...     thetaProductViewer.plot()

The solution will be tested against data that was created with ``steps
= 10`` with a FORTRAN code written by Ryo Kobayashi for phase field
modeling. The following code opens the file `test.gz` extracts the
data and compares it with the `theta` variable.

    >>> import os
    >>> testFile = 'test.gz'
    >>> import examples.phase.impingement.mesh20x20
    >>> import gzip
    >>> filepath = os.path.join(examples.phase.impingement.mesh20x20.__path__[0], testFile)
    >>> filestream = gzip.open(filepath,'r')
    >>> import cPickle
    >>> testData = cPickle.load(filestream)
    >>> filestream.close()
    >>> import Numeric
    >>> if len(testData.flat) == len(Numeric.array(theta)):
    ...     testData = Numeric.reshape(testData, Numeric.array(theta).shape)
    ... else:
    ...     testData = Numeric.resize(testData, Numeric.array(theta).shape)

We iterate the solution in time, plotting as we go if running interactively,

    >>> for i in range(steps):
    ...     theta.updateOld()
    ...     phase.updateOld()
    ...     thetaEq.solve(theta, dt = timeStepDuration)
    ...     phaseEq.solve(phase, dt = timeStepDuration)
    ...     if __name__ == '__main__':
    ...         phaseViewer.plot()
    ... 	thetaProductViewer.plot()
    
The solution is compared against Ryo Kobayashi's test data

    >>> print theta.allclose(testData)
    1

The following code shows how to restart a simulation from some saved
data. First, reset the variables to their original values.

    >>> phase.setValue(0)
    >>> theta.setValue(-pi + 0.0001)
    >>> for a, b, thetaValue in ((0., 0.,  2. * pi / 3.), 
    ...                          (Lx, 0., -2. * pi / 3.), 
    ...                          (0., Lx, -2. * pi / 3. + 0.3), 
    ...                          (Lx, Lx,  2. * pi / 3.)):
    ...     cells = mesh.getCells(filter = cornerCircle, a = a, b = b)
    ...     phase.setValue(1., cells)
    ...     theta.setValue(thetaValue, cells)

Step through half the time steps.

    >>> for i in range(steps / 2):
    ...     theta.updateOld()
    ...     phase.updateOld()
    ...     thetaEq.solve(theta, dt = timeStepDuration)
    ...     phaseEq.solve(phase, dt = timeStepDuration)

We confirm that the solution has not yet converged to that given by 
Ryo Kobayashi's FORTRAN code:

    >>> print theta.allclose(testData)
    0

We save the variables to disk.

    >>> import fipy.tools.dump as dump
    >>> import tempfile
    >>> (f, fileName) = tempfile.mkstemp('.gz')
    >>> dump.write({'phase' : phase, 'theta' : theta, 'thetaEq' : thetaEq, 'phaseEq' : phaseEq}, fileName)
    
and then recall them to test the data pickling mechanism

    >>> data = dump.read(fileName)
    >>> newPhase = data['phase']
    >>> newTheta = data['theta']
    >>> newThetaEq = data['thetaEq']
    >>> newPhaseEq = data['phaseEq']

We clean up the temporary dump file

    >>> import os
    >>> os.close(f)
    >>> os.remove(fileName)

and finish the iterations,

    >>> for i in range(steps = steps / 2):
    ...     newTheta.updateOld()
    ...     newPhase.updateOld()
    ...     newThetaEq.solve(newTheta, dt = timeStepDuration)
    ...     newPhaseEq.solve(newPhase, dt = timeStepDuration)

The solution is compared against Ryo Kobayashi's test data

    >>> print newTheta.allclose(testData)
    1
    
"""
__docformat__ = 'restructuredtext'

def script():
    """
    Return the documentation for this module as a script that can be
    invoked to initialize other scripts.
    """
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus._getScript(__name__)


