#!/usr/bin/env python

## 
 # ###################################################################
 #  PyFiVol - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 3/9/04 {12:02:05 PM} 
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

from fivol.profiler.profiler import Profiler
from fivol.profiler.profiler import calibrate_profiler

from fivol.meshes.grid2D import Grid2D
from fivol.viewers.grid2DGistViewer import Grid2DGistViewer
from fivol.variables.cellVariable import CellVariable

nx = 20
dx = 0.01
L = nx * dx

mesh = Grid2D(
    dx = dx,
    dy = dx,
    nx = nx,
    ny = nx)

var = CellVariable(
    name = "test",
    mesh = mesh)

bottomLeft = mesh.getCells(lambda cell: cell.getCenter()[0] < L/2. and cell.getCenter()[1] < L/2.)
bottomRight = mesh.getCells(lambda cell: cell.getCenter()[0] >= L/2. and cell.getCenter()[1] < L/2.)
topLeft = mesh.getCells(lambda cell: cell.getCenter()[0] < L/2. and cell.getCenter()[1] >= L/2.)
topRight = mesh.getCells(lambda cell: cell.getCenter()[0] >= L/2. and cell.getCenter()[1] >= L/2.)



var.setValue(mesh.getCellCenters()[:,0] * mesh.getCellCenters()[:,1])

if __name__ == '__main__':
    viewer = Grid2DGistViewer(var = var, minVal = 0, maxVal = L * L)
    viewer.plot()
	
    raw_input()
    

##     def bottomRight(cell, cornerCell):
## 	return L - cell.getCenter()[0] == cornerCell.getCenter()[1] and cell.getCenter()[1] == cornerCell.getCenter()[1]
	
##     print mesh.getCellCenters()
    
##     pt = (0.045, 0.045)
##     
##     d = mesh.getPointToCellDistances(pt)
##     
##     print d
## 
##     print Numeric.argmin(d)
## 
##     i = Numeric.argsort(d)
##     
##     print i
##     
##     print (var[i[0]] * var[i[1]] * (d[i[0]] + d[i[1]])) / (var[i[0]] * d[i[0]] + var[i[1]] * d[i[1]])
##     
    print var.getValue(cells = bottomLeft)
##     
##     print bottomRight
    
    fudge = calibrate_profiler(10000)
    profile = Profiler('profile', fudge=fudge)

##     var.setValue(var.getValue(points = [(L - cell.getCenter()[0], cell.getCenter()[1]) for cell in bottomRight]), bottomRight)
##     var.setValue(var.getValue(points = [(cell.getCenter()[0], L - cell.getCenter()[1]) for cell in topLeft]), topLeft)
##     var.setValue(var.getValue(points = [(L - cell.getCenter()[0], L - cell.getCenter()[1]) for cell in topRight]), topRight)

    [var.setValue(var((L - cell.getCenter()[0], cell.getCenter()[1])), [cell]) for cell in bottomRight]
    [var.setValue(var((cell.getCenter()[0], L - cell.getCenter()[1])), [cell]) for cell in topLeft]
    [var.setValue(var((L - cell.getCenter()[0], L - cell.getCenter()[1])), [cell]) for cell in topRight]
    
    profile.stop()
    
##     for cell in bottomRight:
## 	var.setValue(var((L - cell.getCenter()[0], cell.getCenter()[1])), [cell]) 

##     print cornerCells
##     print [mesh.getCells(lambda cell: L - cell.getCenter()[0] == cornerCell.getCenter()[0] and cell.getCenter()[1] == cornerCell.getCenter()[1]) for cornerCell in cornerCells]
##     [var.setValue(var[cornerCell.getID()], mesh.getCells(lambda cell: L - cell.getCenter()[0] == cornerCell.getCenter()[1] and cell.getCenter()[1] == cornerCell.getCenter()[1])) for cornerCell in cornerCells]
    viewer.plot()
	
    raw_input()
    

