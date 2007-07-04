#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "mesh1D.py"
 #                                    created: 11/10/03 {3:23:47 PM}
 #                                last update: 7/3/07 {6:00:11 PM} 
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
 #  2003-11-10 JEG 1.0 original
 # ###################################################################
 ##

r"""

In this example a phase equation is solved in 1 dimension with a
missorientation present. The phase equation is given by:

.. raw:: latex

    $$ \tau_{\phi} \frac{\partial \phi}{\partial t} 
    = \alpha^2 \nabla^2 \phi + \phi ( 1 - \phi ) m_1 ( \phi , T) 
    - 2 s \phi | \nabla \theta | - \epsilon^2 \phi | \nabla \theta |^2 $$

where

.. raw:: latex

    $$ m_1(\phi, T) = \phi - \frac{1}{2} - T \phi ( 1 - \phi ) $$

The initial conditions are:

.. raw:: latex

    \begin{align*}
    \phi &= 1 \qquad \text{for $0 \le x \le L$} \\
    \theta &= \begin{cases}
    1 & \text{for $0 \le x \le L/2$} \\
    0 & \text{for $L/2 < x \le L$}
    \end{cases} \\
    T &= 1 \qquad \text{for $0 \le x \le L$}
    \end{align*}

and boundary conditions

.. raw:: latex

    $\phi = 1$ for $x = 0$ and $x = L$.

.. Further details of the numerical method for this problem can be found in
   "Extending Phase Field Models of Solidification to Polycrystalline
   Materials", J.A. Warren *et al.*, *Acta Materialia*, **51** (2003)
   6035-6058.  
   
Here the phase equation is solved with an explicit technique.

The solution is allowed to evolve for `steps = 100` time steps.

   >>> for step in range(steps):
   ...     phaseEq.solve(phase, dt = timeStepDuration)

The solution is compared with test data. The test data was created
with a FORTRAN code written by Ryo Kobayashi for phase field
modeling. The following code opens the file `mesh1D.gz` extracts the
data and compares it with the `theta` variable.

   >>> import os
   >>> from fipy.tools import dump
   >>> testData = dump.read(os.path.splitext(__file__)[0] + '.gz')
   >>> print phase.allclose(testData)
   1
   
"""
__docformat__ = 'restructuredtext'

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

from fipy.meshes.grid1D import Grid1D
mesh = Grid1D(dx = dx, nx = nx)

from fipy.variables.cellVariable import CellVariable
phase = CellVariable(name = 'PhaseField', mesh = mesh, value = 1.)

from fipy.variables.modularVariable import ModularVariable
theta = ModularVariable(name = 'Theta', mesh = mesh, value = 1.)
theta.setValue(0., where=mesh.getCellCenters()[0] > L / 2.)

from fipy.terms.implicitSourceTerm import ImplicitSourceTerm
mPhiVar = phase - 0.5 + temperature * phase * (1 - phase)
thetaMag = theta.getOld().getGrad().getMag()
implicitSource = mPhiVar * (phase - (mPhiVar < 0))
implicitSource += (2 * s + epsilon**2 * thetaMag) * thetaMag

from fipy.terms.transientTerm import TransientTerm
from fipy.terms.explicitDiffusionTerm import ExplicitDiffusionTerm
phaseEq = TransientTerm(phaseTransientCoeff) == \
          ExplicitDiffusionTerm(alpha**2) \
          - ImplicitSourceTerm(implicitSource) \
          + (mPhiVar > 0) * mPhiVar * phase

if __name__ == '__main__':
   
   import fipy.viewers
   phaseViewer = fipy.viewers.make(vars = phase)
   phaseViewer.plot()
   for step in range(steps):
      phaseEq.solve(phase, dt = timeStepDuration)
      phaseViewer.plot()
   raw_input('finished')
