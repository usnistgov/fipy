#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 4/2/04 {4:06:53 PM} { 5:14:21 PM}
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

import Numeric

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

class AnisotropySystem:

    def __init__(self, n = 40):
        timeStepDuration = 5e-5
        self.steps = 10

        phaseParameters = {
            'tau'                   : 3e-4,
            'time step duration'    : timeStepDuration,
            'epsilon'               : 0.008,
            's'                     : 0.01,
            'alpha'                 : 0.015,
            'anisotropy'            : 0.02,
            'symmetry'              : 4.,
            'kappa 1'               : 0.9,
            'kappa 2'               : 20.
            }

        temperatureParameters = {
            'time step duration'    : timeStepDuration,
            'temperature diffusion' : 2.25,
            'latent heat'           : 1.,
            'heat capacity'         : 1.
            }


        Length = n * 2.5 / 100.
        nx = n
        ny = n
        dx = Length / nx
        dy = Length / ny


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

        self.phaseViewer = Grid2DGistViewer(var = phase)
        self.temperatureViewer = Grid2DGistViewer(var = temperature, minVal = -0.5, maxVal =0.5)
        
        phaseFields = {
            'theta' : theta,
            'temperature' : temperature
            }
        
        temperatureFields = {
            'phase' : phase
            }
        
        def circleCells(cell,L = Length):
            x = cell.getCenter()
            r = L / 4.
            c = (L / 2., L / 2.)
            if (x[0] - c[0])**2 + (x[1] - c[1])**2 < r**2:
                return 1
            else:
                return 0
            
        interiorCells = mesh.getCells(circleCells)
            
        phase.setValue(1.,interiorCells)
        
        phaseEq = PhaseEquation(
            phase,
            mPhi = Type2MPhiVariable,
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
        
        temperatureEq = TemperatureEquation(
            temperature,
            solver = LinearPCGSolver(
            tolerance = 1.e-15, 
            steps = 1000
            ),
            boundaryConditions=(
            FixedFlux(mesh.getExteriorFaces(), 0.),
            ),
            parameters = temperatureParameters,
            fields = temperatureFields
            )

        self.it = Iterator((phaseEq, temperatureEq))

        self.parameters = {
            'it' : self.it,
            'var' : phase,
            'steps' : self.steps
            }

    def getParameters(self):
        return self.parameters

    def run(self):
        self.phaseViewer.plot()
        self.temperatureViewer.plot()

        for i in range(self.steps):
            self.it.timestep(1)
            self.phaseViewer.plot()
            self.temperatureViewer.plot()

if __name__ == '__main__':
    system = AnisotropySystem()
    system.run()
    raw_input()

