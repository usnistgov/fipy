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
from fivol.viewers.grid2DGistViewer import Grid2DGistViewer
from fivol.variables.cellVariable import CellVariable
from fivol.examples.levelSet.levelSetEquation import LevelSetEquation

class LevelSetSystem:

    def __init__(self, nx = 10, ny = 10):

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
            x = cell.getCenter()
            if Lx / 3. < x[0] < 2. * Lx / 3. and Ly / 3. < x[1] < 2. * Ly / 3.:
                return 1
            else:
                return 0

        interiorCells = mesh.getCells(funcInteriorSquare)

        levelSetVariable.setValue(1.,interiorCells)

        self.levelSetEq = LevelSetEquation(levelSetVariable)
            

    def run(self):
        self.levelSetViewer.plot()
        zeroCells = self.levelSetEq.solve()
        for cell in zeroCells:
            print "i: ",cell.getId()
        raw_input()

if __name__ == '__main__':
    system = LevelSetSystem(nx = 5, ny = 5)
    system.run()



