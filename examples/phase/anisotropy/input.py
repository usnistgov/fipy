#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 11/1/04 {11:01:44 AM} { 5:14:21 PM}
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

In this example we solve a coupled phase and temperature equation to model 
solidification, and eventually dendritic growth, from a circular seed in a 2D mesh:
    
>>> numberOfCells = 40
>>> Length = numberOfCells * 2.5 / 100.
>>> nx = numberOfCells
>>> ny = numberOfCells
>>> dx = Length / nx
>>> dy = Length / ny
>>> radius = Length / 4.
>>> from fipy.meshes.grid2D import Grid2D
>>> mesh = Grid2D(dx,dy,nx,ny)
    
Dendritic growth will not be observed with this small test system. If
you wish to see dendritic growth reset the following parameters:
``numberOfCells = 200``, ``steps = 10000``, ``radius = Length / 80``.

The governing equation for the phase field is given by:

.. raw:: latex

    $$ \tau_{\phi} \frac{\partial \phi}{\partial t} 
    = \alpha^2 \nabla^2 \phi + \phi ( 1 - \phi ) m_2 ( \phi , T) 
    - 2 s \phi | \nabla \theta | - \epsilon^2 \phi | \nabla \theta |^2 $$

where

.. raw:: latex

    $$ m_2(\phi, T) 
    = \phi - \frac{1}{2} - \frac{ \kappa_1 }{ \pi } \arctan \left( \kappa_2 T \right) $$
    
and the governing equation for temperature is given by:

.. raw:: latex

    $$ \frac{\partial T}{\partial t} = D_T \nabla^2 T + \frac{\partial \phi}{\partial t} $$

..  Further details of the numerical method for this problem can be found in
    "Extending Phase Field Models of Solidification to Polycrystalline
    Materials", J.A. Warren *et al.*, *Acta Materialia*, **51** (2003) 6035-6058.

Here the phase and temperature equations are solved with an explicit
and implicit technique, respectively.

The parameters for these equations are 

    >>> timeStepDuration = 5e-5
    >>> phaseParameters = {
    ...     'tau'                   : 3e-4,
    ...     'epsilon'               : 0.008,
    ...     's'                     : 0.01,
    ...     'alpha'                 : 0.015,
    ...     'anisotropy'            : 0.02,
    ...     'symmetry'              : 4.,
    ...     'kappa 1'               : 0.9,
    ...     'kappa 2'               : 20.
    ...     }
    >>> temperatureParameters = {
    ...     'timeStepDuration' : timeStepDuration,
    ...     'temperature diffusion' : 2.25,
    ...     'latent heat'           : 1.,
    ...     'heat capacity'         : 1.
    ...     }

The variable `theta` represents the orientation of the crystal. In
this example, it is constant and thus does not affect the solution.

    >>> from fipy.models.phase.theta.modularVariable import ModularVariable
    >>> theta = ModularVariable(
    ...     name = 'Theta',
    ...     mesh = mesh
    ...     )

The `phase` variable is `0` for a liquid and `1` for a solid.  Here we
build an example `phase` variable, initialized as a liquid,

    >>> from fipy.variables.cellVariable import CellVariable
    >>> phase = CellVariable(
    ...     name = 'PhaseField',
    ...     mesh = mesh,
    ...     value = 0.,
    ...     hasOld = 1)

The `hasOld` flag keeps the old value of the variable. This is
necessary for a transient solution. In this example we wish to set up
an interior region that is solid. A value of `1` is assigned to the
`phase` variable on a patch defined by the method:

    >>> def circleCells(cell,L = Length):
    ...     x = cell.getCenter()
    ...     r = radius
    ...     c = (Length / 2., Length / 2.)
    ...     if (x[0] - c[0])**2 + (x[1] - c[1])**2 < r**2:
    ...         return 1
    ...     else:
    ...         return 0
   
This method is passed to `mesh.getCells(filter = circleCells)` which
filters out the required cells.
   
    >>> interiorCells = mesh.getCells(filter = circleCells)           
    >>> phase.setValue(1.,interiorCells)

The temperature field is initialized to a value of `-0.4` throughout:

    >>> temperature = CellVariable(
    ...     name = 'Theta',
    ...     mesh = mesh,
    ...     value = -0.4,
    ...     hasOld = 1
    ...     )
	
The `phase` equation requires a `mPhi` instantiator to represent

.. raw:: latex

   $m_2(\phi, T)$

above

    >>> from fipy.models.phase.phase.type2MPhiVariable import Type2MPhiVariable

The `phase` equation is solved with an iterative conjugate gradient
solver

and requires access to the `theta` and `temperature` variables

    >>> from fipy.models.phase.phase.phaseEquation import buildPhaseEquation
    >>> phaseEq = buildPhaseEquation(
    ...     mPhi = Type2MPhiVariable,
    ...     parameters = phaseParameters,
    ...     theta = theta,
    ...     temperature = temperature,
    ...     phase = phase)
	
The `temperature` equation is also solved with an iterative conjugate
gradient solver and requires access to the `phase` variable

    >>> from fipy.models.phase.temperature.temperatureEquation import buildTemperatureEquation
    >>> temperatureEq = buildTemperatureEquation(
    ...     parameters = temperatureParameters,
    ...     phase = phase)

If we are running this example interactively, we create viewers for
the phase and temperature fields

    >>> if __name__ == '__main__':
    ...     from fipy.viewers.grid2DGistViewer import Grid2DGistViewer
    ...     phaseViewer = Grid2DGistViewer(var = phase)
    ...     temperatureViewer = Grid2DGistViewer(var = temperature, 
    ...                                          minVal = -0.5, maxVal =0.5)
    ...     phaseViewer.plot()
    ...     temperatureViewer.plot()

we iterate the solution in time, plotting as we go if running interactively,

    >>> steps = 10
    >>> for i in range(steps):
    ...     phase.updateOld()
    ...     temperature.updateOld()
    ...     phaseEq.solve(phase, dt = timeStepDuration)
    ...     temperatureEq.solve(temperature, dt = timeStepDuration)
    ...     if i%10 == 0 and __name__ == '__main__':
    ...         phaseViewer.plot()
    ...         temperatureViewer.plot()

The solution is compared with test data. The test data was created for
``steps = 10`` with a FORTRAN code written by Ryo Kobayashi for phase
field modeling. The following code opens the file ``test.gz`` extracts
the data and compares it with the `phase` variable.

   >>> import os
   >>> testFile = 'test.gz'
   >>> import examples.phase.anisotropy
   >>> gzfile = 'gunzip --fast -c < %s/%s'
   >>> gzfile = gzfile%(examples.phase.anisotropy.__path__[0], testFile)
   >>> filestream=os.popen(gzfile,'r')
   >>> import cPickle
   >>> testData = cPickle.load(filestream)
   >>> filestream.close()
   >>> import Numeric
   >>> phase =  Numeric.array(phase)
   >>> testData = Numeric.reshape(testData, phase.shape)
   >>> Numeric.allclose(phase, testData, rtol = 1e-10, atol = 1e-10)
   1
   
"""
__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus.getScript())

    raw_input('finished')

