#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 8/31/04 {1:21:05 PM}
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

In this example we solve the same set of equations as in
`examples/phase/impingemnet/mesh20x20/input.py` but a restart method
is demonstrated. First do `steps / 2` time steps.

   >>> for i in range(steps / 2):        
   ...     it.timestep(dt = timeStepDuration)

Save the variables and recall them to test the mechanism

   >>> import fipy.tools.dump as dump
   >>> dump.write({'phase' : phase, 'theta' : theta, 'mesh' : mesh}, 'data')
   >>> data = dump.read('data')
   >>> newPhase = data['phase']
   >>> newTheta = data['theta']
   >>> newMesh = data['mesh']
   
Rebuild the equations:

   >>> phaseFields = {
   ...     'theta' : newTheta,
   ...     'temperature' : temperature
   ... }
        
   >>> thetaFields = {
   ...     'phase' : newPhase
   ...  }

   >>> newThetaEq = ThetaEquation(
   ...     var = newTheta,
   ...     solver = LinearPCGSolver(
   ...         tolerance = 1.e-15, 
   ...         steps = 2000
   ...     ),
   ...     boundaryConditions = (
   ...         FixedFlux(newMesh.getExteriorFaces(), 0.),
   ...     ),
   ...     parameters = thetaParameters,
   ...     fields = thetaFields
   ... )
           
   >>> newPhaseEq = PhaseEquation(
   ...     var = newPhase,
   ...     mPhi = Type1MPhiVariable,
   ...     solver = LinearPCGSolver(
   ...         tolerance = 1.e-15, 
   ...         steps = 1000
   ...     ),
   ...     boundaryConditions = (
   ...         FixedFlux(newMesh.getExteriorFaces(), 0.),
   ...      ),
   ...      parameters = phaseParameters,
   ...      fields = phaseFields
   ... )

and the iterator:

   >>> newIt = Iterator((newThetaEq, newPhaseEq))
   
Do some more iterations:

   >>> for i in range(5):
   ...     newIt.timestep(dt = timeStepDuration)

and check the results:

   >>> import os
   >>> testFile = 'test.gz'
   >>> import examples.phase.impingement.mesh20x20
   >>> filestream=os.popen('gunzip --fast -c < %s/%s'%(examples.phase.impingement.mesh20x20.__path__[0], testFile),'r')
   >>> import cPickle
   >>> testData = cPickle.load(filestream)
   >>> filestream.close()
   >>> import Numeric
   >>> newTheta =  Numeric.array(newTheta)
   >>> testData = Numeric.reshape(testData, newTheta.shape)
   >>> Numeric.allclose(newTheta, testData, rtol = 1e-10, atol = 1e-10)
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

nx = 20
ny = 20
steps = 10
temperature = 10.

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

pi = Numeric.pi
        
phase = CellVariable(
    name = 'PhaseField',
    mesh = mesh,
    value = 0.,
    hasOld = 1
    )

theta = ModularVariable(
    name = 'Theta',
    mesh = mesh,
    value = -pi + 0.0001,
    hasOld = 1
    )

circleCenters = ((0., 0.), (Lx, 0.), (0., Ly), (Lx, Ly))
thetaValues = (2. * pi / 3., -2. * pi / 3., -2. * pi / 3. + 0.3, 2. * pi / 3.)

for i in range(len(thetaValues)):

    (a, b) = circleCenters[i]
    
    def circle(cell):
        x = cell.getCenter()[0]
        y = cell.getCenter()[1]
        if ((x - a)**2 + (y - b)**2) < (Lx / 2.)**2:
            return 1
        else:
            return 0

    cells = mesh.getCells(filter = circle)

    phase.setValue(1., cells)
    theta.setValue(thetaValues[i], cells)
        

thetaProd = -pi + phase * (theta + pi)
        
phaseViewer = Grid2DGistViewer(var = phase, palette = 'rainbow.gp', minVal = 0., maxVal = 1., grid = 0)
thetaProductViewer = Grid2DGistViewer(var = thetaProd , palette = 'rainbow.gp', minVal = -pi, maxVal = pi, grid = 0)
        
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
    
    phaseViewer.plot()
    thetaProductViewer.plot()

    for i in range(steps):
        it.timestep(dt = timeStepDuration)
        phaseViewer.plot()
        thetaProductViewer.plot()

    raw_input('finished')
