#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "distanceFunctionEquation.py"
 #                                    created: 11/12/03 {10:39:23 AM} 
 #                                last update: 6/3/04 {2:53:54 PM} 
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
import MA

from fipy.equations.equation import Equation
import fipy.tools.vector as vector

class DistanceFunctionEquation(Equation):

    def __init__(self, var):
        
        self.mesh = var.getMesh()
	
	Equation.__init__(
            self,
            var = var,
            terms = (),
            solver = None)

        self.array = self.var.getNumericValue()
        self.initialArray = self.array.copy()
        self.cells = self.mesh.getCells()
        self.resetCells()
        
    def getVar(self):
        self.var.setValue(self.array)
        return self.var

    def solve(self):
        self.resetCells()
        self.setInitialEvaluatedCells()
        trialCellIDs = self.getInitialTrialCells()
        self.iterateOverRemainingCells(trialCellIDs)
        self.var.setValue(self.array)

    def resetCells(self):
        self.array = self.initialArray.copy()
        self.evaluatedFlag = Numeric.zeros(len(self.cells))
        self.trialFlag = Numeric.zeros(len(self.cells))

    def iterateOverRemainingCells(self, trialCellIDs):
        while len(trialCellIDs) > 0:
            cellID = self.getMinAbsCellID(trialCellIDs)
            cell = self.cells[cellID]
            self.evaluatedFlag[cell.getID()] = 1
            trialCellIDs.remove(cellID)
            for adjCellID in cell.getCellToCellIDs():
                if type(adjCellID) is type(1): 
                    if not self.evaluatedFlag[adjCellID]:
                        self.array[adjCellID] = self.getCellValue(self.cells[adjCellID])
                        if not self.trialFlag[adjCellID]:
                            trialCellIDs.append(adjCellID)
                            self.trialFlag[adjCellID] = 1
        
    def setInitialEvaluatedCells(self):
        tmpArray = self.array.copy()

        for cell in self.cells:
            possibleValues = ()

            for index in range(len(cell.getCellToCellIDs())):
                adjCellID = cell.getCellToCellIDs()[index]
                if type(adjCellID) is type(1):
                    val = self.array[cell.getID()]
                    adjVal = self.array[adjCellID]

                    if val * adjVal < 0.0:
                        dAP = cell.getCellToCellDistances()[index]
                        
                        possibleValues += (val * dAP / abs(adjVal - val),)

                        self.evaluatedFlag[cell.getID()] = 1
            
            if self.evaluatedFlag[cell.getID()] == 1:
                if self.array[cell.getID()] > 0.:
                    tmpArray[cell.getID()] = min(possibleValues)
                else:
                    tmpArray[cell.getID()] = max(possibleValues)

        self.array = tmpArray

    def getInitialTrialCells(self):
        trialCellIDs = []
        for cell in self.cells:
            if self.evaluatedFlag[cell.getID()]:
                for adjCellID in cell.getCellToCellIDs():
                    if type(adjCellID) is type(1):
                        if not self.evaluatedFlag[adjCellID]:
                            if not self.trialFlag[adjCellID]:
                                self.trialFlag[adjCellID] = 1
                                trialCellIDs.append(adjCellID)
                                self.array[adjCellID] = self.getCellValue(self.cells[adjCellID])

        return trialCellIDs

    def getAdjacentEvaluatedIndices(self, cell):
        adjIndices = []
        for index in range(len(cell.getCellToCellIDs())):
            id = cell.getCellToCellIDs()[index]
            if type(id) is type(1):
                if self.evaluatedFlag[id]:
                    adjIndices.append(index)
        return adjIndices

    def getCellValue(self, cell):
        adjIndices = self.getAdjacentEvaluatedIndices(cell)

        if len(adjIndices) == 0:
            raise Exception
        elif len(adjIndices) == 1:
            return self.evaluateOneCell(cell, adjIndices)
        elif len(adjIndices) == 2:
            return self.evaluateTwoCells(cell, adjIndices)
        elif len(adjIndices) == 3:
            return self.evaluateThreeCells(cell, adjIndices)
        elif len(adjIndices) == 4:
            return self.evaluateFourCells(cell, adjIndices)
        else:
            raise Exception
        
    def getMinCellID(self, cellIDs):
        minCellID = cellIDs[0]
        for cellID in cellIDs[1:]:
            if self.array[cellID] < self.array[minCellID]:
                minCellID = cellID

        return minCellID

    def getMaxCellID(self, cellIDs):
        maxCellID = cellIDs[0]
        for cellID in cellIDs[1:]:
            if self.array[cellID] > self.array[maxCellID]:
                maxCellID = cellID

        return maxCellID

    def getMinAbsCellID(self, cellIDs):
        numCellIDs = Numeric.array(cellIDs)
        arr = Numeric.take(self.array, numCellIDs)

        return numCellIDs[Numeric.argsort(Numeric.absolute(arr))[0]]

##        minCellID = cells[0]
##        for cell in cells[1:]:
##            if abs(self.array[cell.getID()]) < abs(self.array[minCell.getID()]):
##                minCell = cell

##        return minCell
            
    def evaluateOneCell(self, cell, indices):
        dAP = cell.getCellToCellDistances()[indices[0]]
        phiP = self.array[cell.getID()]
        phiA = self.array[cell.getCellToCellIDs()[indices[0]]]
        if phiP > 0.0:
            return phiA + dAP
        else:
            return phiA - dAP

    def evaluateTwoCells(self, cell, indices):
        n0 = cell.getNormal(indices[0])
        n1 = cell.getNormal(indices[1])
        if vector.sqrtDot(n0, n1) > 0.99:
            vals = (self.evaluateOneCell(cell, (indices[0],)), self.evaluateOneCell(cell, (indices[1],)))
            if self.array[cell.getID()] > 0.:
                return min(vals)
            else:
                return max(vals)
        else:
            return self._evaluateTwoCells(cell, indices[0], indices[1])

    def evaluateThreeCells(self, cell, indices):
        if self.array[cell.getID()] < 0.:
            cellIDs = MA.take(cell.getCellToCellIDs(), indices)
            minCellID = self.getMinCellID(cellIDs)
            filledList = list(MA.filled(cell.getCellToCellIDs(), value = -1))
            indices.remove(filledList.index(minCellID))
            return self.evaluateTwoCells(cell, indices)
        else:
            cellIDs = MA.take(cell.getCellToCellIDs(), indices)
            maxCellID = self.getMaxCellID(cellIDs)
            filledList = list(MA.filled(cell.getCellToCellIDs(), value = -1))
            indices.remove(filledList.index(maxCellID))
            return self.evaluateTwoCells(cell, indices)

    def evaluateFourCells(self, cell, indices):

        if self.array[cell.getID()] < 0.:
            cellIDs = Numeric.take(cell.getCellToCellIDs(), indices)
            minCellID = self.getMinCellID(cellIDs)
            indices.remove(list(cell.getCellToCellIDs()).index(minCellID))
            return self.evaluateThreeCells(cell, indices)
        else:
            cellIDs = Numeric.take(cell.getCellToCellIDs(), indices)
            maxCellID = self.getMaxCellID(cellIDs)
            indices.remove(list(cell.getCellToCellIDs()).index(maxCellID))
            return self.evaluateThreeCells(cell, indices)

    def _evaluateTwoCells(self, cell, index1, index2):
        d1 = cell.getCellToCellDistances()[index1]
        d2 = cell.getCellToCellDistances()[index2]
        phi = self.array[cell.getID()]
        phi1 = self.array[cell.getCellToCellIDs()[index1]]
        phi2 = self.array[cell.getCellToCellIDs()[index2]]
        aa = d1**2 + d2**2
        bb = -2 * (phi1 * d2**2 + phi2 * d1**2)
        cc = phi1**2 * d2**2 + phi2**2 * d1**2 - d1**2 * d2**2
        sqr = Numeric.sqrt(bb**2 - 4 * aa * cc)
        if phi > 0.:
            return (-bb + sqr) / 2. / aa
        else:
            return (-bb - sqr) / 2. / aa
        
