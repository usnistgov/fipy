#!/usr/bin/env python

## 
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "input.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 01/08/04 { 4:14:57 PM}
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

from meshes.grid2D import Grid2D
from examples.phase.phase.type1PhaseEquation import Type1PhaseEquation
from solvers.linearPCGSolver import LinearPCGSolver
from boundaryConditions.fixedValue import FixedValue
from boundaryConditions.fixedFlux import FixedFlux
from iterators.iterator import Iterator
from viewers.grid2DGistViewer import Grid2DGistViewer
from variables.cellVariable import CellVariable
from examples.phase.theta.modularVariable import ModularVariable
from examples.phase.temperature.temperatureEquation import TemperatureEquation
from examples.phase.theta.thetaEquation import ThetaEquation
from __future__ import nested_scopes
import Numeric

class ThetaSystem:

    def __init__(self, nx = 40, ny = 40, initialConditions = None):
        timeStepDuration = 0.02
        timeStepDuration = 1.
        self.steps = 10

        sharedPhaseThetaParameters = {
            'epsilon'               : 0.008,
            's'                     : 0.01,
            'anisotropy'            : 0.0,
            'alpha'                 : 0.015,
            'symmetry'              : 4.
            }

        phaseParameters = {
            'tau'                   : 0.1,
            'time step duration'    : timeStepDuration,
            'kappa 1'               : 0.9,
            'kappa 2'               : 20.
            }

        thetaParameters = {
            'time step duration'    : timeStepDuration,
            'small value'           : 1e-6,
            'beta'                  : 1e5,
            'mu'                    : 1e3,
            'tau'                   : 3e-5,
            'gamma'                 : 1e3 
            }

        for key in sharedPhaseThetaParameters.keys():
            phaseParameters[key] = sharedPhaseThetaParameters[key]
            thetaParameters[key] = sharedPhaseThetaParameters[key]
            
##        Lx = 2.5 * nx / 100.
##        Ly = 2.5 * ny / 100.            
##        dx = Lx / nx
##        dy = Ly / ny
        dx = 1.
        dy = 1.
        Lx = 10.
        Ly = 10.
        mesh = Grid2D(dx,dy,nx,ny)

        phase = CellVariable(
            name = 'PhaseField',
            mesh = mesh,
            value = 0.
            )
        
        theta = ModularVariable(
            name = 'Theta',
            mesh = mesh,
            value = 0.,
            hasOld = 0
            )
        
        self.phaseViewer = Grid2DGistViewer(var = phase)
        self.thetaViewer = Grid2DGistViewer(var = theta)
        
        phaseFields = {
            'theta' : theta,
            'temperature' : -0.4
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

        self.it = Iterator((phaseEq, thetaEq))

        self.parameters = {
            'it' : self.it,
            'phase' : phase,
            'theta' : theta,
            'steps' : self.steps
            }

    def getParameters(self):
        return self.parameters

    def run(self):
        self.phaseViewer.plot()
        self.thetaViewer.plot()
        raw_input()
        for i in range(self.steps):
            self.it.timestep(1)
            self.phaseViewer.plot()
            self.thetaViewer.plot()

if __name__ == '__main__':
    def getRightCells(cell, Lx = 1., Ly = 1.):
        if cell.getCenter()[0] > Lx / 2.:
            return 1

    def getAllCells(cell, Lx = 1., Ly = 1.):
        return 1.
    
    initialConditions = (
        { 'phase value' : 1., 'theta value' : 1., 'func' : getAllCells },
        { 'phase value' : 1., 'theta value' : 0., 'func' : getRightCells }        
        )
    
    system = ThetaSystem(nx = 10, ny = 1,  initialConditions = initialConditions)
    system.run()
    raw_input()

