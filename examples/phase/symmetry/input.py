#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 5/5/04 {6:41:41 PM} 
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

from fipy.tools.profiler.profiler import Profiler
from fipy.tools.profiler.profiler import calibrate_profiler

from fipy.meshes.grid2D import Grid2D
from fipy.viewers.grid2DGistViewer import Grid2DGistViewer
from fipy.variables.cellVariable import CellVariable


class SymmetrySystem:
    def __init__(self, N = 20, L = 1.):
        self.N = N
        nx = self.N
        dx = L / N
        L = N * dx

        mesh = Grid2D(
            dx = dx,
            dy = dx,
            nx = N,
            ny = N)

        self.var = CellVariable(
            name = "test",
            mesh = mesh)

        self.var.setValue(mesh.getCellCenters()[:,0] * mesh.getCellCenters()[:,1])

        def func(cell, Length = None):
            return cell.getCenter()[0] < Length / 2. and cell.getCenter()[1] < Length / 2.

        bottomLeftCells = mesh.getCells(filter = func,  Length = L)
        bottomRightCells = ()
        topLeftCells = ()
        topRightCells = ()
        
        for cell in bottomLeftCells:
            x, y = cell.getCenter()
            bottomRightCells += (mesh.getNearestCell((L - x, y)),)            
            topRightCells += (mesh.getNearestCell((L - x , L - y)),)
            topLeftCells += (mesh.getNearestCell((x , L - y)),)

        self.orderedCells = (bottomRightCells, topRightCells, topLeftCells)
        self.symmetryCells = bottomLeftCells
        
        self.viewer = Grid2DGistViewer(var = self.var, minVal = 0, maxVal = L * L)

    def plot(self):
        self.viewer.plot()
        
    def run(self):
        for cellSet in self.orderedCells:
            for i in range(len(cellSet)):
                id = self.symmetryCells[i].getID()
                idOther = cellSet[i].getID()
                self.var[idOther] = self.var[id]

    def getVar(self):
        return self.var

if __name__ == '__main__':
    system = SymmetrySystem(N = 100)
    system.plot()
    raw_input('press key to continue')
    fudge = calibrate_profiler(10000)
    profile = Profiler('profile.txt', fudge=fudge)
    system.run()
    profile.stop()
    system.plot()
    raw_input('press key to continue')
    



