#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 5/5/04 {6:41:20 PM} { 5:14:21 PM}
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
 # protection and is in the public domain.  PFM is an experimental
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

"""

In this example we solve a coupled phase and temperature equation to
model solidification and eventually dendritic growth. Dendritic growth will
not be observed with this small test system. If you wish to see dendritic growth
reset the following parameters:

   >>> numberOfCells = 200
   >>> steps = 10000
   >>> radius = Length / 80.

The following equations are solved, phase equation:

.. raw:: latex

    $$ \\tau_{\\phi} \\frac{\\partial \\phi}{\\partial t} = \\alpha^2 \\nabla^2 \\phi + \\phi ( 1 - \\phi ) m_2 ( \\phi , T) - 2 s \\phi | \\nabla \\theta | - \\epsilon^2 \\phi | \\nabla \\theta |^2 $$

where

.. raw:: latex

    $$ m_2(\\phi, T) = \\phi - \\frac{1}{2} - \\frac{ \\kappa_1 }{ \\pi } \\arctan \\left( \\kappa_2 T \\right) $$
    
and the temperature equation is given by:

.. raw:: latex

    $$ \\frac{\\partial T}{\\partial t} = D_T \\nabla^2 T + \\frac{\\partial \\phi}{\\partial t} $$

Further details of the numerical method for this problem can be found
in Extending Phase Field Models of Solidification to Polycrystalline
Materials, J.A. Warren et al., Acta Materialia, 51 (2003) 6035-6058.
Here the phase and temperature equations are solved with an explicit
and implicit technique respectively.

The parameters for these equations are given in `phaseParameters` and
`temperatureParameters`. The variable `theta` represents the
orientation of the crystal. Here it is constant and thus does not
affect the solution. The `phase` variable is 0 for a liquid and 1 for
a solid.  Here we build an example `phaseVar` invoking with a zero
value,

   >>> phaseVar = CellVariable(
   ...     name = 'PhaseField',
   ...     mesh = mesh,
   ...     value = 0.,
   ...     hasOld = 1)

The `hasOld` flag keeps the old value of the variable. This is
necessary for a transient solution. In this example we wish to set up
an interior region that is solid. A value of 1 if given to the `phase`
variable on a patch defined by the method `circleCells`. This method
is passed to `mesh.getCells(filter = circleCells)` which filters out
the required cells.

   >>> def circleCells(cell,L = Length):
   ...     x = cell.getCenter()
   ...     r = radius
   ...     c = (Length / 2., Length / 2.)
   ...     if (x[0] - c[0])**2 + (x[1] - c[1])**2 < r**2:
   ...         return 1
   ...     else:
   ...         return 0
   >>> interiorCells = mesh.getCells(filter = circleCells)           
   >>> phaseVar.setValue(1.,interiorCells)

The `phaseEquation` requires a `mPhi` instantiator. Here we use
`Type2MPhiVariable` as outlined in the above equations.` To compare
with the test result the problem is iterated for `steps = 10` time
steps.

   >>> steps = 10
   >>> for i in range(steps):
   ...     it.timestep(dt = timeStepDuration)

The solution is compared with test data. The test data was created
with a FORTRAN code written by Ryo Kobayshi for phase field
modeling. The following code opens the file `test.gz` extracts the
data and compares it with the `phase` variable.

   >>> import os
   >>> testFile = 'test.gz'
   >>> import examples.phase.anisotropy
   >>> filestream=os.popen('gunzip --fast -c < %s/%s'%(examples.phase.anisotropy.__path__[0], testFile),'r')
   >>> import cPickle
   >>> testData = cPickle.load(filestream)
   >>> filestream.close()
   >>> import Numeric
   >>> phase =  Numeric.array(phase)
   >>> testData = Numeric.reshape(testData, phase.shape)
   >>> Numeric.allclose(phase, testData, rtol = 1e-10, atol = 1e-10)
   1
   
"""

from fipy.meshes.grid2D import Grid2D
from fipy.models.phase.phase.type2MPhiVariable import Type2MPhiVariable
from fipy.models.phase.phase.phaseEquation import PhaseEquation
from fipy.solvers.linearPCGSolver import LinearPCGSolver
from fipy.boundaryConditions.fixedValue import FixedValue
from fipy.boundaryConditions.fixedFlux import FixedFlux
from fipy.iterators.iterator import Iterator
from fipy.viewers.grid2DGistViewer import Grid2DGistViewer
from fipy.variables.cellVariable import CellVariable
from fipy.models.phase.theta.modularVariable import ModularVariable
from fipy.models.phase.temperature.temperatureEquation import TemperatureEquation

timeStepDuration = 5e-5
steps = 10

phaseParameters = {
    'tau'                   : 3e-4,
    'epsilon'               : 0.008,
    's'                     : 0.01,
    'alpha'                 : 0.015,
    'anisotropy'            : 0.02,
    'symmetry'              : 4.,
    'kappa 1'               : 0.9,
    'kappa 2'               : 20.
    }

temperatureParameters = {
    'timeStepDuration' : timeStepDuration,
    'temperature diffusion' : 2.25,
    'latent heat'           : 1.,
    'heat capacity'         : 1.
    }


numberOfCells = 40
Length = numberOfCells * 2.5 / 100.
nx = numberOfCells
ny = numberOfCells
dx = Length / nx
dy = Length / ny
radius = Length / 4.

mesh = Grid2D(dx,dy,nx,ny)

phase = CellVariable(
    name = 'PhaseField',
    mesh = mesh,
    value = 0.,
    hasOld = 1
    )
        
theta = ModularVariable(
    name = 'Theta',
    mesh = mesh,
    value = 0.,
    hasOld = 0
    )
        
temperature = CellVariable(
    name = 'Theta',
    mesh = mesh,
    value = -0.4,
    hasOld = 1
    )

phaseViewer = Grid2DGistViewer(var = phase)
temperatureViewer = Grid2DGistViewer(var = temperature, minVal = -0.5, maxVal =0.5)
        
phaseFields = {
    'theta' : theta,
    'temperature' : temperature
    }
        
temperatureFields = {
    'phase' : phase
    }
        
def circleCells(cell):
    x = cell.getCenter()
    r = radius
    c = (Length / 2., Length / 2.)
    if (x[0] - c[0])**2 + (x[1] - c[1])**2 < r**2:
        return 1
    else:
        return 0
            
interiorCells = mesh.getCells(filter = circleCells)
            
phase.setValue(1.,interiorCells)
        
phaseEq = PhaseEquation(
    phase,
    mPhi = Type2MPhiVariable,
    solver = LinearPCGSolver(
    tolerance = 1.e-15,
    steps = 1000
    ),
    boundaryConditions=(FixedFlux(mesh.getExteriorFaces(), 0.),),
    parameters = phaseParameters,
    fields = phaseFields
    )
        
temperatureEq = TemperatureEquation(
    temperature,
    solver = LinearPCGSolver(
    tolerance = 1.e-15, 
    steps = 1000
    ),
    boundaryConditions=(FixedFlux(mesh.getExteriorFaces(), 0.),),
    parameters = temperatureParameters,
    fields = temperatureFields
    )

it = Iterator((phaseEq, temperatureEq))

if __name__ == '__main__':
    
    phaseViewer.plot()
    temperatureViewer.plot()

    for i in range(steps):
        it.timestep(dt = timeStepDuration)
        if i%10 == 0:
            phaseViewer.plot()
            temperatureViewer.plot()

    raw_input('finished')

