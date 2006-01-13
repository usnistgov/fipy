#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 1/12/06 {9:29:15 PM}
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
 #  2003-11-17 JEG 1.0 original
 # ###################################################################
 ##

r"""

In this example we solve a coupled phase and orientation equation on a
one dimensional grid.

    >>> nx = 40
    >>> Lx = 2.5 * nx / 100.
    >>> dx = Lx / nx
    >>> from fipy.meshes.grid1D import Grid1D
    >>> mesh = Grid1D(dx, nx)
	
This problem simulates the wet boundary that forms between grains of different 
orientations. The phase equation is given by

.. raw:: latex

    $$ \tau_{\phi} \frac{\partial \phi}{\partial t} 
    = \alpha^2 \nabla^2 \phi + \phi ( 1 - \phi ) m_1 ( \phi , T) 
    - 2 s \phi | \nabla \theta | - \epsilon^2 \phi | \nabla \theta |^2 $$

where

.. raw:: latex

    $$ m_1(\phi, T) = \phi - \frac{1}{2} - T \phi ( 1 - \phi ) $$

and the orientation equation is given by

.. raw:: latex

    $$ P(\epsilon | \nabla \theta |) \tau_{\theta} \phi^2 
    \frac{\partial \theta}{\partial t} 
    = \nabla \cdot \left[ \phi^2 \left( \frac{s}{| \nabla \theta |} 
    + \epsilon^2 \right) \nabla \theta \right] $$

where

.. raw:: latex

    $$ P(w) = 1 - \exp{(-\beta w)} + \frac{\mu}{\epsilon} \exp{(-\beta w)} $$

The initial conditions for this problem are set such that

.. raw:: latex

   $\phi = 1$ for $0 \le x \le L_x$ and 
   
   $$
   \theta = \begin{cases}
   1 & \text{for $0 \le x < L_x / 2$,} \\
   0 & \text{for $L_x / 2 \le x \le L_x$.}
   \end{cases}
   $$

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

    >>> from fipy.variables.cellVariable import CellVariable
    >>> phase = CellVariable(
    ...     name = 'PhaseField',
    ...     mesh = mesh,
    ...     value = 1.
    ...     )

Because `theta`

.. raw:: latex

    is an $S^1$-valued variable (i.e. it maps to the circle) and thus
    intrinsically has $2\pi$-peridocity,

we must use `ModularVariable` instead of a `CellVariable`. A
`ModularVariable` confines `theta` to

.. raw:: latex

    $-\pi < \theta \le \pi$ by adding or subtracting $2\pi$ where
    necessary and by defining a new

subtraction operator between two angles.

    >>> from fipy.variables.modularVariable import ModularVariable
    >>> theta = ModularVariable(
    ...     name = 'Theta',
    ...     mesh = mesh,
    ...     value = 1.,
    ...     hasOld = 1
    ...     )

The left and right halves of the domain are given different orientations.
    
    >>> theta.setValue(0., where=mesh.getCellCenters()[...,0] > Lx / 2.)

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

The `phase` equation is constructed.

    >>> phaseEq = TransientTerm(phaseTransientCoeff) == ExplicitDiffusionTerm(alpha**2) \
    ...                                                 - ImplicitSourceTerm(implicitSource) \
    ...                                                 + (mPhiVar > 0) * mPhiVar * phase

The `theta` equation is built in the following way. The details for
this equation are fairly involved, see J.A. Warren *et al.*. The main
detail is that a source must be added to correct for the
discretization of `theta` on the circle.

    >>> phaseMod = phase + ( phase < thetaSmallValue ) * thetaSmallValue
    >>> phaseModSq = phaseMod * phaseMod
    >>> expo = epsilon * beta * theta.getGrad().getMag()
    >>> expo = (expo < 100.) * (expo - 100.) + 100.
    >>> from fipy.tools import numerix
    >>> pFunc = 1. + numerix.exp(-expo) * (mu / epsilon - 1.)

    >>> phaseFace = phase.getArithmeticFaceValue()
    >>> phaseSq = phaseFace * phaseFace
    >>> gradMag = theta.getFaceGrad().getMag()
    >>> eps = 1. / gamma / 10.
    >>> gradMag += (gradMag < eps) * eps
    >>> IGamma = (gradMag > 1. / gamma) * (1 / gradMag - gamma) + gamma
    >>> diffusionCoeff = phaseSq * (s * IGamma + epsilon**2)

The source term requires the evaluation of the face gradient without
the modular operator. A method of `ModularVariable`, `getFaceGradNoMod()`,
evelautes the gradient without modular arithmetic.

    >>> thetaGradDiff = theta.getFaceGrad() - theta.getFaceGradNoMod()
    >>> sourceCoeff = (diffusionCoeff * thetaGradDiff).getDivergence()

Finally the `theta` equation can be constructed.

    >>> from fipy.terms.implicitDiffusionTerm import ImplicitDiffusionTerm
    >>> thetaEq = TransientTerm(thetaTransientCoeff * phaseModSq * pFunc) == \
    ...           ImplicitDiffusionTerm(diffusionCoeff) \
    ...           + sourceCoeff

If the example is run interactively, we create viewers for the phase
and orientation variables.

    >>> if __name__ == '__main__':
    ...     from fipy import viewers
    ...     phaseViewer = viewers.make(vars = phase, 
    ...                                limits = {'datamin': 0., 'datamax': 1.})
    ...     from fipy.tools.numerix import pi
    ...     thetaProductViewer = viewers.make(vars = theta,
    ...                                       limits = {'datamin': -pi, 
    ...                                                 'datamax': pi})
    ...     phaseViewer.plot()
    ...     thetaProductViewer.plot()

we iterate the solution in time, plotting as we go if running interactively,

    >>> steps = 10
    >>> for i in range(steps):
    ...     theta.updateOld()
    ...     phase.updateOld()
    ...     thetaEq.solve(theta, dt = timeStepDuration)
    ...     phaseEq.solve(phase, dt = timeStepDuration)
    ...     if __name__ == '__main__':
    ...         phaseViewer.plot()
    ... 	thetaProductViewer.plot()

The solution is compared with test data. The test data was created
with ``steps = 10`` with a FORTRAN code written by Ryo Kobayashi for
phase field modeling. The following code opens the file `test.gz`
extracts the data and compares it with the `theta` variable.

   >>> import os
   >>> testFile = 'test.gz'
   >>> import examples.phase.impingement.mesh40x1
   >>> filepath = os.path.join(examples.phase.impingement.mesh40x1.__path__[0],
   ...                         testFile)
   >>> import fipy.tools.dump as dump
   >>> testData = dump.read(filepath)
   >>> print theta.allclose(testData)
   1
"""
__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())

    raw_input('finished')
