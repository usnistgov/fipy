#!/usr/bin/env python

## 
 # ###################################################################
 #  PyFiVol - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 01/27/04 { 1:23:41 PM}
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

class LevelSetSystem:

    def __init__(self):

        nx = 10
        ny = 10
        dx = 1.
        dy = 1.

        mesh = Grid2D(dx,dy,nx,ny)

        levelSetVariable = CellVariable(
            name = 'level set variable',
            mesh = mesh,
            value = -1.
            )

        self.levelSetViewer = Grid2DGistViewer(var = levelSetVariable, palette = 'rainbow.gp', minVal = -1., maxVal = 1.)

        Lx = nx * dx
        Ly = ny * dy

        def funcInteriorSquare(cell, Lx = Lx, Ly = Ly):
            cell.getCenter()
        

        interiorCells = mesh .getCells(funcInteriorSquare

        levelSetVariable.setValue(1.,interiorCells)


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
        self.phaseViewer.plot()
        self.thetaViewer.plot()
##        raw_input()
        for i in range(self.steps):
            self.it.timestep(1)
            self.phaseViewer.plot()
            self.thetaViewer.plot()
##            raw_input()



