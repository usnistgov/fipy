#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 11/1/04 {11:02:02 AM}
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

In this example we solve a coupled phase and orientation equation on a one 
dimensional grid

    >>> nx = 40
    >>> ny = 1
    >>> Lx = 2.5 * nx / 100.
    >>> Ly = 2.5 * ny / 100.            
    >>> dx = Lx / nx
    >>> dy = Ly / ny
    >>> from fipy.meshes.grid2D import Grid2D
    >>> mesh = Grid2D(dx,dy,nx,ny)
	
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
    >>> phaseParameters = {
    ...    'tau'                   : 0.1,
    ...     }
    >>> thetaParameters = {
    ...     'small value'           : 1e-6,
    ...     'beta'                  : 1e5,
    ...     'mu'                    : 1e3,
    ...     'tau'                   : 0.01,
    ...     'gamma'                 : 1e3 
    ...     }

with the shared parameters

    >>> sharedPhaseThetaParameters = {
    ...     'epsilon'               : 0.008,
    ...     's'                     : 0.01,
    ...     'anisotropy'            : 0.0,
    ...     'alpha'                 : 0.015,
    ...     'symmetry'              : 4.
    ...     }
    >>> for key in sharedPhaseThetaParameters.keys():
    ...     phaseParameters[key] = sharedPhaseThetaParameters[key]
    ...     thetaParameters[key] = sharedPhaseThetaParameters[key]

The system is held isothermal at

    >>> temperature = 1.

and is initially solid everywhere

    >>> from fipy.variables.cellVariable import CellVariable
    >>> phase = CellVariable(
    ...     name = 'PhaseField',
    ...     mesh = mesh,
    ...     value = 1.
    ...     )

The left and right halves of the domain are given different orientations

    >>> from fipy.models.phase.theta.modularVariable import ModularVariable
    >>> theta = ModularVariable(
    ...     name = 'Theta',
    ...     mesh = mesh,
    ...     value = 1.,
    ...     hasOld = 1
    ...     )
    >>> theta.setValue(0., mesh.getCells(filter = lambda cell: cell.getCenter()[0] > Lx / 2.))

The `phase` equation requires a `mPhi` instantiator to represent

.. raw:: latex

   $m_1(\phi, T)$

above

    >>> from fipy.models.phase.phase.type1MPhiVariable import Type1MPhiVariable

and requires access to the `theta` and `temperature` variables

    >>> from fipy.models.phase.phase.phaseEquation import buildPhaseEquation
    >>> phaseEq = buildPhaseEquation(
    ...     phase = phase,
    ...     mPhi = Type1MPhiVariable,
    ...     parameters = phaseParameters,
    ...     theta = theta,
    ...     temperature = temperature)

The `theta` equation is also solved with an iterative conjugate gradient solver  
and requires access to the `phase` variable

    >>> from fipy.models.phase.theta.thetaEquation import buildThetaEquation
    >>> thetaEq = buildThetaEquation(
    ...     theta = theta,
    ...     parameters = thetaParameters,
    ...     phase = phase)

If the example is run interactively, we create viewers for the phase and 
orientation variables. Rather than viewing the raw orientation, which is not 
meaningful in the liquid phase, we weight the orientation by the phase

    >>> if __name__ == '__main__':
    ...     from fipy.viewers.grid2DGistViewer import Grid2DGistViewer
    ...     phaseViewer = Grid2DGistViewer(var = phase, palette = 'rainbow.gp', 
    ...                                    minVal = 0., maxVal = 1., grid = 0)
    ...     from Numeric import pi
    ...     thetaProd = -pi + phase * (theta + pi)
    ...     thetaProductViewer = Grid2DGistViewer(var = thetaProd , palette = 'rainbow.gp', 
    ...                                           minVal = -pi, maxVal = pi, grid = 0)
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

The solution is compared with test data. The test data was created with 
``steps = 10`` with a FORTRAN code written by Ryo Kobayashi for phase field
modeling. The following code opens the file `test.gz` extracts the
data and compares it with the `theta` variable.

   >>> import os
   >>> testFile = 'test.gz'
   >>> import examples.phase.impingement.mesh40x1
   >>> import gzip
   >>> filepath = os.path.join(examples.phase.impingement.mesh40x1.__path__[0], testFile)
   >>> filestream = gzip.open(filepath,'r')
   >>> import cPickle
   >>> testData = cPickle.load(filestream)
   >>> filestream.close()
   >>> import Numeric
   >>> testData = Numeric.reshape(testData, Numeric.array(theta).shape)
   >>> theta.allclose(testData)
   1
"""
__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus.getScript())

    raw_input('finished')
