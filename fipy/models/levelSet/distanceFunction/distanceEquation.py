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

class DistanceEquation(Equation):

    def __init__(self, var):
        
        self.mesh = var.getMesh()
	
	Equation.__init__(
            self,
            var = var,
            terms = (),
            solver = None)

    def solve(self):
        setValueFlag = self._calcInterfaceValues()
        setValueFlag = self._calcInitialTrialValues(setValueFlag)
##        self._calcRemainingValues(setValueFlag)

    def _calcInterfaceValues(self):
        """

        Sets the values in cells at the interface (cells that have a neighbour
        of the opposite sign) to the shortest signed distance.

           >>> from fipy.meshes.grid2D import Grid2D
           >>> mesh = Grid2D(dx = 1., dy = 1., nx = 2, ny = 1)
           >>> from distanceVariable import DistanceVariable
           >>> var = DistanceVariable(mesh = mesh, value = (-1, 1))
           >>> DistanceEquation(var)._calcInterfaceValues()
           [1,1,]
           >>> print var
           [-0.5, 0.5,]

           >>> mesh = Grid2D(dx = 1., dy = 1., nx = 3, ny = 3)
           >>> from distanceVariable import DistanceVariable
           >>> var = DistanceVariable(mesh = mesh, value =
           ...     (-1, -1, -1, -1, 1, -1, -1, -1, -1))
           >>> DistanceEquation(var)._calcInterfaceValues()
           [0,1,0,1,1,1,0,1,0,]
           >>> answer = Numeric.array((-1,-0.5,-1,-0.5,0.5,-0.5,-1,-0.5,-1))
           >>> Numeric.allclose(answer, Numeric.array(var))
           1

           >>> mesh = Grid2D(dx = .5, dy = 2., nx = 3, ny = 3)
           >>> from distanceVariable import DistanceVariable
           >>> var = DistanceVariable(mesh = mesh, value =
           ...     (-1, -1, -1, -1, 1, -1, -1, -1, -1))
           >>> DistanceEquation(var)._calcInterfaceValues()
           [0,1,0,1,1,1,0,1,0,]
           >>> answer = Numeric.array((-1,-1,-1,-0.25,0.25,-0.25,-1,-1,-1))
           >>> Numeric.allclose(answer, Numeric.array(var))
           1

        """
        
        N = self.mesh.getNumberOfCells()
        M = self.mesh.getMaxFacesPerCell()
        
        cellZeroFlag = Numeric.take(self.var.getInterfaceFlag(), self.mesh.getCellFaceIDs())
        
        dAP = self.mesh.getCellToCellDistances()
        cellToCellIDs = self.mesh.getCellToCellIDsFilled()
        phiAdj = Numeric.take(self.var, cellToCellIDs)
        phi = Numeric.resize(Numeric.repeat(self.var, M),(N, M))


        distance = MA.masked_values(phi * dAP / abs(phi - phiAdj) * cellZeroFlag, 0)

        argmins = MA.argmin(abs(distance), axis = 1)

        argmins = Numeric.transpose(Numeric.array((Numeric.arange(N),argmins)))

        argmins = argmins[:,0] * M + argmins[:,1]

        distance = MA.take(distance.flat, argmins)
        
        cellFlag = Numeric.sum(cellZeroFlag, axis = 1)

        self.var.setValue(MA.where(cellFlag > 0, distance, self.var))
        return Numeric.where(cellFlag > 0, 1, 0)

    def _calcQuadratic(self, phi1, phi2, n1, n2, d1, d2):
        """

        Calculates the distance function from two neighbouring values.
        It should work on any given mesh.

           >>> from fipy.meshes.grid2D import Grid2D
           >>> mesh = Grid2D(dx = 1., dy = 1., nx = 1, ny = 1)
           >>> from distanceVariable import DistanceVariable
           >>> var = DistanceVariable(mesh = mesh, value = 1)
           >>> eqn = DistanceEquation(var)
           >>> sqrt = Numeric.sqrt(2)
           >>> n1 = Numeric.array((1,0))
           >>> n2 = Numeric.array((0,1))
           >>> answer = Numeric.array((1 + 1/ sqrt, 1 - 1 / sqrt))
           >>> Numeric.allclose(answer, eqn._calcQuadratic(1, 1, n1, n2, 1, 1))
           1
           >>> n1 = Numeric.array((-1,-1)) / sqrt
           >>> n2 = Numeric.array((1,-1)) / sqrt
           >>> d = 1
           >>> answer = Numeric.array((1 + d, 1 - d))
           >>> Numeric.allclose(answer, eqn._calcQuadratic(1, 1, n1, n2, d * sqrt, d * sqrt))
           1
           >>> n1 = Numeric.array((0,-1))
           >>> n2 = Numeric.array((-1,-1)) / sqrt
           >>> d = 0.3
           >>> answer = Numeric.array((1 + d, 1 - d))
           >>> Numeric.allclose(answer, eqn._calcQuadratic(1, 1, n1, n2, d, d * sqrt))
           1
           
        """

        dotProd = d1 * d2 * Numeric.dot(n1, n2)
        crossProd = d1 * d2 * (n1[0] * n2[1] - n1[1] * n2[0])
        dsq = d1**2 + d2**2 - 2 * dotProd

        top = -phi1 * (dotProd - d2**2) - phi2 * (dotProd - d1**2)
        sqrt = Numeric.sqrt(crossProd**2 *(dsq - (phi1 - phi2)**2))
        
        return Numeric.array(((top + sqrt) / dsq, (top - sqrt) / dsq))
        
    def _calcInitialTrialValues(self, setValueFlag):
        """

        Sets that variable value at one trial value.

           >>> from fipy.meshes.grid2D import Grid2D
           >>> mesh = Grid2D(dx = 1., dy = 1., nx = 1, ny = 3)
           >>> from distanceVariable import DistanceVariable
           >>> var = DistanceVariable(mesh = mesh, value = (-1, 1, 1))
           >>> eqn = DistanceEquation(var)
           >>> setValueFlag = eqn._calcInterfaceValues()
           >>> setValueFlag = eqn._calcInitialTrialValues(setValueFlag)
           >>> Numeric.allclose(Numeric.array((-.5, .5, 1.5)), Numeric.array(var))
           1

           >>> mesh = Grid2D(dx = 1., dy = 1., nx = 2, ny = 2)
           >>> var = DistanceVariable(mesh = mesh, value = (-1, 1, 1, 1))
           >>> eqn = DistanceEquation(var)
           >>> setValueFlag = eqn._calcInterfaceValues()
           >>> eqn._calcInitialTrialValues(setValueFlag)
           [1,1,1,2,]
           >>> sqrt = Numeric.sqrt(2)
           >>> Numeric.allclose(Numeric.array((-.5, .5, .5, .5 + 1 / sqrt)), Numeric.array(var))
           1

           >>> mesh = Grid2D(dx = .5, dy = 2., nx = 3, ny = 3)
           >>> var = DistanceVariable(mesh = mesh, value = (-1, -1, -1, -1, 1, -1, -1, -1, -1))
           >>> eqn = DistanceEquation(var)
           >>> setValueFlag = eqn._calcInterfaceValues()
           >>> eqn._calcInitialTrialValues(setValueFlag)
           [2,1,2,1,1,1,2,1,2,]
           >>> v = -1.40771446
           >>> Numeric.allclose(Numeric.array((v, -1, v, -.25, .25, -.25, v, -1, v)), Numeric.array(var), atol = 1e-6)
           1
        """
        N = self.mesh.getNumberOfCells()
        M = self.mesh.getMaxFacesPerCell()
        cellToCellIDs = self.mesh.getCellToCellIDsFilled()

        adjacentCellFlag = Numeric.take(setValueFlag, cellToCellIDs)

        sumAdjacentCellFlag = Numeric.sum(adjacentCellFlag, axis = 1)
        
        setValueFlag = Numeric.where(sumAdjacentCellFlag == 0,
                                     setValueFlag,
                                     Numeric.where(setValueFlag == 1,
                                                   setValueFlag,
                                                   2))

        

        for cellID in range(N):
            if setValueFlag[cellID] == 2:
                self.var[cellID] = self._calcTrialValue(cellID, setValueFlag)

        return setValueFlag
     
    def _calcTrialValue(self, cellID, setValueFlag):

        cellToCellIDs = self.mesh.getCellToCellIDsFilled()[cellID]
        dAP = self.mesh.getCellToCellDistances()[cellID]
        phiAdj = Numeric.take(self.var, cellToCellIDs)
        values = MA.array(phiAdj, mask = (Numeric.take(setValueFlag, cellToCellIDs)!=1))

        argsort = MA.argsort(abs(values))

        values = MA.take(values, argsort)
        dAP = Numeric.take(self.mesh.getCellToCellDistances()[cellID], argsort)
        normals = Numeric.take(self.mesh.getCellNormals()[cellID], argsort)

        sign = -1
        if self.var[cellID] > 0:
            sign = 1


        NSetValues = len(values.mask()) - Numeric.sum(values.mask())

        if NSetValues == 0:
            raise Error
        elif NSetValues == 1:
            return values[0] + sign * dAP[0]
        elif NSetValues == 2:
            quad = self._calcQuadratic(values[0], values[1], normals[0], normals[1], dAP[0], dAP[1])
            if sign > 0:
                return quad[0]
            else:
                return quad[1]

def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__": 
    _test() 



##    def resetCells(self):
##        self.array = self.initialArray.copy()
##        self.evaluatedFlag = Numeric.zeros(len(self.cells))
##        self.trialFlag = Numeric.zeros(len(self.cells))

##    def iterateOverRemainingCells(self, trialCellIDs):
##        while len(trialCellIDs) > 0:
##            cellID = self.getMinAbsCellID(trialCellIDs)
##            cell = self.cells[cellID]
##            self.evaluatedFlag[cell.getID()] = 1
##            trialCellIDs.remove(cellID)
##            for adjCellID in cell.getCellToCellIDs():
##                if type(adjCellID) is type(1): 
##                    if not self.evaluatedFlag[adjCellID]:
##                        self.array[adjCellID] = self.getCellValue(self.cells[adjCellID])
##                        if not self.trialFlag[adjCellID]:
##                            trialCellIDs.append(adjCellID)
##                            self.trialFlag[adjCellID] = 1
        
##    def setInitialEvaluatedCells(self):
##        tmpArray = self.array.copy()

##        for cell in self.cells:
##            possibleValues = ()

##            for index in range(len(cell.getCellToCellIDs())):
##                adjCellID = cell.getCellToCellIDs()[index]
##                if type(adjCellID) is type(1):
##                    val = self.array[cell.getID()]
##                    adjVal = self.array[adjCellID]

##                    if val * adjVal < 0.0:
##                        dAP = cell.getCellToCellDistances()[index]
                        
##                        possibleValues += (val * dAP / abs(adjVal - val),)

##                        self.evaluatedFlag[cell.getID()] = 1
            
##            if self.evaluatedFlag[cell.getID()] == 1:
##                if self.array[cell.getID()] > 0.:
##                    tmpArray[cell.getID()] = min(possibleValues)
##                else:
##                    tmpArray[cell.getID()] = max(possibleValues)

##        self.array = tmpArray

##    def getInitialTrialCells(self):
##        trialCellIDs = []
##        for cell in self.cells:
##            if self.evaluatedFlag[cell.getID()]:
##                for adjCellID in cell.getCellToCellIDs():
##                    if type(adjCellID) is type(1):
##                        if not self.evaluatedFlag[adjCellID]:
##                            if not self.trialFlag[adjCellID]:
##                                self.trialFlag[adjCellID] = 1
##                                trialCellIDs.append(adjCellID)
##                                self.array[adjCellID] = self.getCellValue(self.cells[adjCellID])

##        return trialCellIDs

##    def getAdjacentEvaluatedIndices(self, cell):
##        adjIndices = []
##        for index in range(len(cell.getCellToCellIDs())):
##            id = cell.getCellToCellIDs()[index]
##            if type(id) is type(1):
##                if self.evaluatedFlag[id]:
##                    adjIndices.append(index)
##        return adjIndices

##    def getCellValue(self, cell):
##        adjIndices = self.getAdjacentEvaluatedIndices(cell)

##        if len(adjIndices) == 0:
##            raise Exception
##        elif len(adjIndices) == 1:
##            return self.evaluateOneCell(cell, adjIndices)
##        elif len(adjIndices) == 2:
##            return self.evaluateTwoCells(cell, adjIndices)
##        elif len(adjIndices) == 3:
##            return self.evaluateThreeCells(cell, adjIndices)
##        elif len(adjIndices) == 4:
##            return self.evaluateFourCells(cell, adjIndices)
##        else:
##            raise Exception
        
##    def getMinCellID(self, cellIDs):
##        minCellID = cellIDs[0]
##        for cellID in cellIDs[1:]:
##            if self.array[cellID] < self.array[minCellID]:
##                minCellID = cellID

##        return minCellID

##    def getMaxCellID(self, cellIDs):
##        maxCellID = cellIDs[0]
##        for cellID in cellIDs[1:]:
##            if self.array[cellID] > self.array[maxCellID]:
##                maxCellID = cellID

##        return maxCellID

##    def getMinAbsCellID(self, cellIDs):
##        numCellIDs = Numeric.array(cellIDs)
##        arr = Numeric.take(self.array, numCellIDs)

##        return numCellIDs[Numeric.argsort(Numeric.absolute(arr))[0]]

####        minCellID = cells[0]
####        for cell in cells[1:]:
####            if abs(self.array[cell.getID()]) < abs(self.array[minCell.getID()]):
####                minCell = cell

####        return minCell
            
##    def evaluateOneCell(self, cell, indices):
##        dAP = cell.getCellToCellDistances()[indices[0]]
##        phiP = self.array[cell.getID()]
##        phiA = self.array[cell.getCellToCellIDs()[indices[0]]]
##        if phiP > 0.0:
##            return phiA + dAP
##        else:
##            return phiA - dAP

##    def evaluateTwoCells(self, cell, indices):
##        n0 = cell.getNormal(indices[0])
##        n1 = cell.getNormal(indices[1])
##        if vector.sqrtDot(n0, n1) > 0.99:
##            vals = (self.evaluateOneCell(cell, (indices[0],)), self.evaluateOneCell(cell, (indices[1],)))
##            if self.array[cell.getID()] > 0.:
##                return min(vals)
##            else:
##                return max(vals)
##        else:
##            return self._evaluateTwoCells(cell, indices[0], indices[1])

##    def evaluateThreeCells(self, cell, indices):
##        if self.array[cell.getID()] < 0.:
##            cellIDs = MA.take(cell.getCellToCellIDs(), indices)
##            minCellID = self.getMinCellID(cellIDs)
##            filledList = list(MA.filled(cell.getCellToCellIDs(), value = -1))
##            indices.remove(filledList.index(minCellID))
##            return self.evaluateTwoCells(cell, indices)
##        else:
##            cellIDs = MA.take(cell.getCellToCellIDs(), indices)
##            maxCellID = self.getMaxCellID(cellIDs)
##            filledList = list(MA.filled(cell.getCellToCellIDs(), value = -1))
##            indices.remove(filledList.index(maxCellID))
##            return self.evaluateTwoCells(cell, indices)

##    def evaluateFourCells(self, cell, indices):

##        if self.array[cell.getID()] < 0.:
##            cellIDs = Numeric.take(cell.getCellToCellIDs(), indices)
##            minCellID = self.getMinCellID(cellIDs)
##            indices.remove(list(cell.getCellToCellIDs()).index(minCellID))
##            return self.evaluateThreeCells(cell, indices)
##        else:
##            cellIDs = Numeric.take(cell.getCellToCellIDs(), indices)
##            maxCellID = self.getMaxCellID(cellIDs)
##            indices.remove(list(cell.getCellToCellIDs()).index(maxCellID))
##            return self.evaluateThreeCells(cell, indices)

##    def _evaluateTwoCells(self, cell, index1, index2):
##        d1 = cell.getCellToCellDistances()[index1]
##        d2 = cell.getCellToCellDistances()[index2]
##        phi = self.array[cell.getID()]
##        phi1 = self.array[cell.getCellToCellIDs()[index1]]
##        phi2 = self.array[cell.getCellToCellIDs()[index2]]
##        aa = d1**2 + d2**2
##        bb = -2 * (phi1 * d2**2 + phi2 * d1**2)
##        cc = phi1**2 * d2**2 + phi2**2 * d1**2 - d1**2 * d2**2
##        sqr = Numeric.sqrt(bb**2 - 4 * aa * cc)
##        if phi > 0.:
##            return (-bb + sqr) / 2. / aa
##        else:
##            return (-bb - sqr) / 2. / aa
        

