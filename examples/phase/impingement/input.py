#!/usr/bin/env python

## 
 # ###################################################################
 #  PyFiVol - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 1/16/04 {12:00:06 PM}
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

from __future__ import nested_scopes

import Numeric

from fivol.meshes.grid2D import Grid2D
from fivol.examples.phase.phase.type1PhaseEquation import Type1PhaseEquation
from fivol.solvers.linearPCGSolver import LinearPCGSolver
from fivol.boundaryConditions.fixedValue import FixedValue
from fivol.boundaryConditions.fixedFlux import FixedFlux
from fivol.iterators.iterator import Iterator
from fivol.viewers.grid2DGistViewer import Grid2DGistViewer
from fivol.variables.cellVariable import CellVariable
from fivol.examples.phase.theta.modularVariable import ModularVariable
from fivol.examples.phase.temperature.temperatureEquation import TemperatureEquation
from fivol.examples.phase.theta.thetaEquation import ThetaEquation

class ImpingementSystem:

    def __init__(self, nx = 100, ny = 100, initialConditions = None, steps = 10, drivingForce = 1.):
        timeStepDuration = 0.02
        self.steps = steps

        sharedPhaseThetaParameters = {
            'epsilon'               : 0.008,
            's'                     : 0.01,
            'anisotropy'            : 0.0,
            'alpha'                 : 0.015,
            'symmetry'              : 4.
            }

        phaseParameters = {
            'tau'                   : 0.1,
            'time step duration'    : timeStepDuration
            }

        thetaParameters = {
            'time step duration'    : timeStepDuration,
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

##        Lx = 10
##        Ly = 10
##        nx = 10
##        ny = 10
##        dx = 1.
##        dy = 1.

        mesh = Grid2D(dx,dy,nx,ny)

        phase = CellVariable(
            name = 'PhaseField',
            mesh = mesh,
            value = 0.
            )
        
        theta = ModularVariable(
            name = 'Theta',
            mesh = mesh,
            value = 0.
            )

        pi = Numeric.pi
        self.phaseViewer = Grid2DGistViewer(var = phase, palette = 'rainbow.gp', minVal = 0., maxVal = 1.)
        self.thetaViewer = Grid2DGistViewer(var = theta, palette = 'rainbow.gp', minVal = -pi, maxVal = pi)
        
        phaseFields = {
            'theta' : theta,
            'temperature' : drivingForce
            }
        
        thetaFields = {
            'phase' : phase
            }

        for initialCondition in initialConditions:
            func = initialCondition['func']
            def funcIn(cell, Lx = Lx, Ly = Ly):
                return func(cell, Lx = Lx, Ly = Ly)
            cells = mesh.getCells(funcIn)
            phase.setValue(initialCondition['phase value'],cells)
            theta.setValue(initialCondition['theta value'],cells)
            
        thetaEq = ThetaEquation(
            var = theta,
            solver = LinearPCGSolver(
            tolerance = 1.e-15, 
            steps = 1000
            ),
            boundaryConditions=(
            FixedFlux(mesh.getExteriorFaces(), 0.),
            ),
            parameters = thetaParameters,
            fields = thetaFields
            )
        
        phaseEq = Type1PhaseEquation(
            var = phase,
            solver = LinearPCGSolver(
            tolerance = 1.e-15, 
            steps = 1000
            ),
            boundaryConditions=(
            FixedFlux(mesh.getExteriorFaces(), 0.),
            ),
            parameters = phaseParameters,
            fields = phaseFields
            )
        
        self.it = Iterator((thetaEq, phaseEq))

        self.parameters = {
            'it' : self.it,
            'phase' : phase,
            'theta' : theta,
            'steps' : self.steps
            }

    def getParameters(self):
        return self.parameters

    def run(self):
        self.phaseViewer.plot(fileName = 'phase.ps')
        raw_input('written file, press key to continue')
        self.thetaViewer.plot()

        for i in range(self.steps):
            self.it.timestep(1)
            self.phaseViewer.plot()
            self.thetaViewer.plot()

        self.thetaViewer.plot(fileName = 'theta.ps')


