#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 10/15/04 {11:46:16 AM}
 #  Author: Jonathan Guyer
 #  E-mail: guyer@nist.gov
 #  Author: Daniel Wheeler
 #  E-mail: daniel.wheeler@nist.gov
 #    mail: NIST
 #     www: http://ctcms.nist.gov
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

In this example we solve a coupled phase and orientation equation. It
is a 1D problem that simulates the wet boundary that forms between
grains of different orientations. The phase equation is given by:

.. raw:: latex

    $$ \tau_{\phi} \frac{\partial \phi}{\partial t} 
    = \alpha^2 \nabla^2 \phi + \phi ( 1 - \phi ) m_1 ( \phi , T) 
    - 2 s \phi | \nabla \theta | - \epsilon^2 \phi | \nabla \theta |^2 $$

where

.. raw:: latex

    $$ m_1(\phi, T) = \phi - \frac{1}{2} - T \phi ( 1 - \phi ) $$

an the orientation equation is given by,

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
    0 & \text{for $L_x / 2 \ge x \le L_x$.}
    \end{cases}
    $$

.. Further details of the numerical method for this problem can be found in
   "Extending Phase Field Models of Solidification to Polycrystalline
   Materials", J.A. Warren *et al.*, *Acta Materialia*, **51** (2003)
   6035-6058.  

Here the phase and orientation equations are solved with an
explicit and implicit technique respectively.

The `phaseEquation` requires an `mPhi` instantiator. Here we use
`Type1MPhiVariable` as outlined in the above equations. To compare
with the test result the problem is iterated for `steps = 10` time
steps.

   >>> for i in range(steps):        
   ...     it.timestep(dt = timeStepDuration)

The solution is compared with test data. The test data was created
with a FORTRAN code written by Ryo Kobayshi for phase field
modeling. The following code opens the file `test.gz` extracts the
data and compares it with the `theta` variable.

   >>> import os
   >>> testFile = 'test.gz'
   >>> import examples.phase.impingement.mesh40x1
   >>> gzfile = 'gunzip --fast -c < %s/%s'
   >>> gzfile = gzfile%(examples.phase.impingement.mesh40x1.__path__[0], testFile)
   >>> filestream=os.popen(gzfile,'r')
   >>> import cPickle
   >>> testData = cPickle.load(filestream)
   >>> filestream.close()
   >>> import Numeric
   >>> theta =  Numeric.array(theta)
   >>> testData = Numeric.reshape(testData, theta.shape)
   >>> Numeric.allclose(theta, testData, rtol = 1e-10, atol = 1e-10)
   1


"""
__docformat__ = 'restructuredtext'

##from __future__ import nested_scopes

import Numeric

from fipy.meshes.grid2D import Grid2D
from fipy.models.phase.phase.type1MPhiVariable import Type1MPhiVariable
from fipy.models.phase.phase.phaseEquation import PhaseEquation
from fipy.solvers.linearPCGSolver import LinearPCGSolver
from fipy.boundaryConditions.fixedValue import FixedValue
from fipy.boundaryConditions.fixedFlux import FixedFlux
from fipy.iterators.iterator import Iterator
from fipy.viewers.grid2DGistViewer import Grid2DGistViewer
from fipy.variables.cellVariable import CellVariable
from fipy.models.phase.theta.modularVariable import ModularVariable
from fipy.models.phase.theta.thetaEquation import ThetaEquation


nx = 40
ny = 1
steps = 10
temperature = 1.

timeStepDuration = 0.02

sharedPhaseThetaParameters = {
    'epsilon'               : 0.008,
    's'                     : 0.01,
    'anisotropy'            : 0.0,
    'alpha'                 : 0.015,
    'symmetry'              : 4.
    }

phaseParameters = {
    'tau'                   : 0.1,
    }

thetaParameters = {
    'small value'           : 1e-6,
    'beta'                  : 1e5,
    'mu'                    : 1e3,
    'tau'                   : 0.01,
    'gamma'                 : 1e3 
    }

for key in sharedPhaseThetaParameters.keys():
    phaseParameters[key] = sharedPhaseThetaParameters[key]
    thetaParameters[key] = sharedPhaseThetaParameters[key]
            
Lx = 2.5 * nx / 100.
Ly = 2.5 * ny / 100.            
dx = Lx / nx
dy = Ly / ny

mesh = Grid2D(dx,dy,nx,ny)
        
phase = CellVariable(
    name = 'PhaseField',
    mesh = mesh,
    value = 1.
    )

theta = ModularVariable(
    name = 'Theta',
    mesh = mesh,
    value = 1.,
    hasOld = 1
    )

def getRightCells(cell):
    if cell.getCenter()[0] > Lx / 2.:
        return 1

theta.setValue(0., mesh.getCells(filter = getRightCells))

pi = Numeric.pi
thetaProd = -pi + phase * (theta + pi)
        
phaseFields = {
    'theta' : theta,
    'temperature' : temperature
    }
        
thetaFields = {
    'phase' : phase
    }

thetaEq = ThetaEquation(
    var = theta,
    solver = LinearPCGSolver(
    tolerance = 1.e-15, 
    steps = 2000
    ),
    boundaryConditions = (
    FixedFlux(mesh.getExteriorFaces(), 0.),
    ),
    parameters = thetaParameters,
    fields = thetaFields
    )
        
phaseEq = PhaseEquation(
    var = phase,
    mPhi = Type1MPhiVariable,
    solver = LinearPCGSolver(
    tolerance = 1.e-15, 
    steps = 1000
    ),
    boundaryConditions = (
    FixedFlux(mesh.getExteriorFaces(), 0.),
    ),
    parameters = phaseParameters,
    fields = phaseFields
    )
        
it = Iterator((thetaEq, phaseEq))

if __name__ == '__main__':
    phaseViewer = Grid2DGistViewer(var = phase, palette = 'rainbow.gp', minVal = 0., maxVal = 1., grid = 0)
    thetaProductViewer = Grid2DGistViewer(var = thetaProd , palette = 'rainbow.gp', minVal = -pi, maxVal = pi, grid = 0)
    phaseViewer.plot()
    thetaProductViewer.plot()

    for i in range(steps):
        it.timestep(dt = timeStepDuration)
        phaseViewer.plot()
        thetaProductViewer.plot()

    raw_input('finished')
