#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  PyFiVol - Python-based finite volume PDE solver
 # 
 #  FILE: "levelSetEquation.py"
 #                                    created: 11/12/03 {10:39:23 AM} 
 #                                last update: 1/16/04 {11:37:22 AM} 
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
 #  2003-11-12 JEG 1.0 original
 # ###################################################################
 ##

import Numeric

from fivol.equations.equation import Equation

class LevelSetEquation(Equation):
    """
    Level set equation is implicit.
    """    
    def __init__(self, var):
        
        self.mesh = var.getMesh()
	
	Equation.__init__(
            self,
            var = var,
            terms = (),
            solver = None)

    def solve(self):
        cells = self.mesh.getCells()
        zeroCells = self.getZeroCells(cells)
        self.setZeroValues(zeroCells)

        ## find bounding cells to the evaluatedCells
##        boundingCells = self.getBoundingCells(zeroCells)

##        while(len(boundingCells) != 0):

##        ## get MinimumCell
##            minimumCell = self.getMinimumCell(boundingCells)

##        ## set value of minimum cell
##            self.setCellValue(minimumCell)

##        ## add the extra cell
##            boundingCells, zeroCells = self.updateBoundingCells(minimumCell, boundingCells, zeroCells)

    def getZeroCells(self, cells):
        zeroCells = ()
        array = self.var.getNumericValue()
        for cell in cells:
            id = cell.getId()
            zeroCell = ()
            for adjacentCell in cell.getBoundingCells():                
                adjacentId = adjacentCell.getId()
                if array[id] * array[adjacentId]<0.:
                    zeroCell = (cell,)
            zeroCells += zeroCell
        return zeroCells

    def setZeroCellValues(self, zeroCells):
        varOld = var.copy()
        for cell in zeroCells:
            minCell1 = cell.getMinimumCell(cell.getAdjacentCells())
            minCell2 = cell.getMinimumCell(cell.getAdjacentCells()-minCell1)
            value = varOld.getValue(cell)
            gradx,grady = cell.getGradient(minCell1, minCell2)
            value = value / Numeric.sqrt(gradx * gradx, grady * grady)
            var.setValue(value, cell)

##    def getBoundingCells(self, zeroCells):
##        boundingCells = ()
##        for cell in cells:
##            boundingCell = ()
##            if cell !in zeroCells:
##                for adjacentCell in cell.getAdjacentCells():
##                    if adjacentCell in zeroCells:
##                        boundingCell = (cell,)
##            boundingCells += boundingCell

    def getMinimumCell(cells, var):
        minVal = var(cells[0].getId())
        for cell in cells[1:]:
            minVal = min(minVal, var(cell.getId()))

        return minVal



