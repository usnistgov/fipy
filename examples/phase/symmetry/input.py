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

"""

   >>> testResult = Numeric.zeros((N / 2, N / 2), 'd')
   >>> bottomRight = Numeric.zeros((N / 2, N / 2), 'd')
   >>> topLeft = Numeric.zeros((N / 2, N / 2), 'd')
   >>> topRight = Numeric.zeros((N / 2, N / 2), 'd')
   >>> for j in range(N / 2):
   ...     for i in range(N / 2):
   ...         x = dx * (i + 0.5)
   ...         y = dx * (j + 0.5)
   ...         testResult[i, j] = x * y
   ...         bottomRight[i,j] = var((L - x, y)).getNumericValue()
   ...         topLeft[i,j] = var((x, L - y)).getNumericValue()
   ...         topRight[i,j] = var((L - x, L - y)).getNumericValue()
   >>> Numeric.allclose(testResult, bottomRight, atol = 1e-10)
   1
   >>> Numeric.allclose(testResult,topLeft, atol = 1e-10)
   1
   >>> Numeric.allclose(testResult,topRight, atol = 1e-10)
   1
   
"""

import Numeric

from fipy.meshes.grid2D import Grid2D
from fipy.viewers.grid2DGistViewer import Grid2DGistViewer
from fipy.variables.cellVariable import CellVariable

N = 20
L = 1.
dx = L / N
dy = L / N

mesh = Grid2D(
    dx = dx,
    dy = dy,
    nx = N,
    ny = N)

var = CellVariable(
    name = "test",
    mesh = mesh)

var.setValue(mesh.getCellCenters()[:,0] * mesh.getCellCenters()[:,1])

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
    
    orderedCells = (bottomRightCells, topRightCells, topLeftCells)
    symmetryCells = bottomLeftCells
        
viewer = Grid2DGistViewer(var = var, minVal = 0, maxVal = L * L)

for cellSet in orderedCells:
    for i in range(len(cellSet)):
        id = symmetryCells[i].getID()
        idOther = cellSet[i].getID()
        var[idOther] = var[id]

if __name__ == '__main__':

    viewer.plot()
    raw_input('finished')
    



