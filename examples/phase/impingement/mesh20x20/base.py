#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 10/26/04 {9:00:00 PM} 
 #                                last update: 3/7/05 {1:31:55 PM}
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

    >>> numberOfElements = parse('--numberOfElements', action = 'store', type = 'int', default = 400)
    >>> numberOfSteps = parse('--numberOfSteps', action = 'store', type = 'int', default = 10)

    >>> steps = numberOfSteps
    >>> import Numeric
    >>> nx = int(Numeric.sqrt(numberOfElements))
    >>> ny = nx
    >>> Lx = 2.5 * nx / 100.
    >>> Ly = 2.5 * ny / 100.
    >>> dx = Lx / nx
    >>> dy = Ly / ny
    >>> from fipy.meshes.grid2D import Grid2D
    >>> mesh = Grid2D(dx,dy,nx,ny)

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
    >>> phaseParameters = {
    ...    'tau'                   : 0.1,
    ...    'time step duration'    : timeStepDuration
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

This time, the system is held isothermal at

    >>> temperature = 10.
    
and is initialized to liquid everywhere

    >>> from fipy.variables.cellVariable import CellVariable
    >>> phase = CellVariable(
    ...     name = 'PhaseField',
    ...     mesh = mesh,
    ...     value = 0.
    ...     )

The orientation is initialized to a uniform value to denote the randomly 
oriented liquid phase

    >>> from Numeric import pi
    >>> from fipy.models.phase.theta.modularVariable import ModularVariable
    >>> theta = ModularVariable(
    ...     name = 'Theta',
    ...     mesh = mesh,
    ...     value = -pi + 0.0001,
    ...     hasOld = 1
    ...     )

Four different solid circular domains are created at each corner of the domain 
with appropriate orientations

    >>> def cornerCircle(cell, a = 1., b = 1., Lx = Lx, Ly = Ly):
    ...     x = cell.getCenter()[0]
    ...     y = cell.getCenter()[1]
    ...     if ((x - a)**2 + (y - b)**2) < (Lx / 2.)**2:
    ...         return 1
    ...     else:
    ...         return 0
    
    >>> for a, b, thetaValue in ((0., 0.,  2. * pi / 3.), 
    ...                          (Lx, 0., -2. * pi / 3.), 
    ...                          (0., Ly, -2. * pi / 3. + 0.3), 
    ...                          (Lx, Ly,  2. * pi / 3.)):
    ...     cells = mesh.getCells(filter = cornerCircle, a = a, b = b)
    ...     phase.setValue(1., cells)
    ...     theta.setValue(thetaValue, cells)

The `phase` equation requires a `mPhi` instantiator to represent

.. raw:: latex

   $m_1(\phi, T)$

above

    >>> from fipy.models.phase.phase.type1MPhiVariable import Type1MPhiVariable

and requires access to the `theta` and `temperature` variables

    >>> from fipy.models.phase.phase.phaseEquation import buildPhaseEquation
    >>> phaseEq = buildPhaseEquation(
    ...     mPhi = Type1MPhiVariable,
    ...     phase = phase,
    ...     theta = theta,
    ...     temperature = temperature,
    ...     parameters = phaseParameters)

The `theta` equation is also solved with an iterative conjugate gradient solver  
and requires access to the `phase` variable

    >>> from fipy.models.phase.theta.thetaEquation import buildThetaEquation
    >>> thetaEq = buildThetaEquation(
    ...     theta = theta,
    ...     parameters = thetaParameters,
    ...     phase = phase)
    
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
    
The preceding initialization steps are used in the next few examples.
"""
__docformat__ = 'restructuredtext'

import optparse

class MyOptionParser(optparse.OptionParser):
    def error(self, msg):
        if 'no such option' in msg:
            pass
        else:
            optparse.OptionParser(self, msg)


def script():
    """
    Return the documentation for this module as a script that can be
    invoked to initialize other scripts.
    """
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.getScript(__name__)


