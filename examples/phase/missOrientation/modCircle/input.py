#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 11/10/03 {3:23:47 PM}
 #                                last update: 4/5/05 {8:08:51 PM} 
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

In this example a phase equation is solved in one dimension with a
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
   \phi &= 1 \qquad \forall x  \\
   \theta &= \begin{cases}
   2 \pi / 3 & \text{for $(x - L / 2)^2 + (y - L / 2)^2 > (L / 4)^2$} \\
   -2 \pi / 3 & \text{for $(x - L / 2)^2 + (y - L / 2)^2 \le (L / 4)^2$}
   \end{cases} \\
   T &= 1 \qquad \forall x 
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
   ...     phase.updateOld()
   ...     eq.solve(phase, dt = timeStepDuration)

The solution is compared with test data. The test data was created
with a FORTRAN code written by Ryo Kobayashi for phase field
modeling. The following code opens the file `test.gz` extracts the
data and compares it with the `theta` variable.

   >>> import os
   >>> import examples.phase.missOrientation.modCircle
   >>> import gzip
   >>> filepath = os.path.join(examples.phase.missOrientation.modCircle.__path__[0], 'test.gz')
   >>> filestream = gzip.open(filepath,'r')
   >>> import cPickle
   >>> testData = cPickle.load(filestream)
   >>> filestream.close()
   >>> import Numeric
   >>> testData = Numeric.reshape(testData, Numeric.array(phase).shape)
   >>> print phase.allclose(testData)
   1
   
"""
__docformat__ = 'restructuredtext'

import Numeric

from fipy.meshes.grid2D import Grid2D
from fipy.models.phase.phase.type1MPhiVariable import Type1MPhiVariable
from fipy.models.phase.theta.modularVariable import ModularVariable
from fipy.variables.cellVariable import CellVariable
import fipy.viewers

steps = 100
timeStepDuration = 0.02
L = 1.5
nx = 100
ny = 100

phaseParameters={
   'tau' :        0.1,
   'epsilon' :    0.008,
   's' :          0.01,
   'alpha' :      0.015,
   'c2':          0.0,
   'anisotropy':  0.,
   'symmetry':    4.
   }
      
valueLeft=1.
valueRight=1.
      
dx = L / nx
dy = L / ny

mesh = Grid2D(dx, dy, nx, ny)

phase = CellVariable(
   name = 'PhaseField',
   mesh = mesh,
   value = 1.
   )
      
theta = ModularVariable(
   name = 'Theta',
   mesh = mesh,
   value = 2. * Numeric.pi / 3.,
   hasOld = 0
   )

def circleFunc(cell):
   r = L / 4.
   c = (L / 2., L / 2.)
   x = cell.getCenter()
   return (x[0] - c[0])**2 + (x[1] - c[1])**2 < r**2

circleCells = mesh.getCells(filter = circleFunc)

theta.setValue(-2. * Numeric.pi / 3., circleCells)
      
temperature = 1.

from fipy.models.phase.phase.phaseEquation import buildPhaseEquation
eq = buildPhaseEquation(
   phase = phase,
   mPhi = Type1MPhiVariable,
   temperature = temperature,
   theta = theta,
   parameters = phaseParameters
   )

if __name__ == '__main__':
   phaseViewer = fipy.viewers.make(vars = phase)
   phaseViewer.plot()
   for step in range(steps):
      phase.updateOld()
      eq.solve(phase, dt = timeStepDuration)
      phaseViewer.plot()
   raw_input('finished')
