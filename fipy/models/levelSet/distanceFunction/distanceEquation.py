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

"""

A `DistanceEquation` object solves the equation,

.. raw:: latex

    $$ | \\nabla \\phi | = 1 $$

using the fast marching method with an initial condition defined by
the zero level set.

Currently the solution is first order, This suffices for initial
conditions with straight edges (e.g. trenches in
electrodeposition). The method should work for unstructured 2D grids
but testing on unstructured grids is untested thus far. This is a 2D
implementation as it stands. Extending to 3D should be relatively
simple.

Here we will define a few test cases. Firstly a 1D test case

   >>> from fipy.meshes.grid2D import Grid2D
   >>> mesh = Grid2D(dx = .5, dy = .2, nx = 8, ny = 1)
   >>> from distanceVariable import DistanceVariable
   >>> var = DistanceVariable(mesh = mesh, value = (-1, -1, -1, -1, 1, 1, 1, 1))
   >>> eqn = DistanceEquation(var)
   >>> eqn.solve()
   >>> answer = (-1.75, -1.25, -.75, -0.25, 0.25, 0.75, 1.25, 1.75) 
   >>> Numeric.allclose(answer, Numeric.array(var))
   1

"""
__docformat__ = 'restructuredtext'

import Numeric
import MA

from fipy.equations.equation import Equation
import fipy.tools.vector as vector

class DistanceEquation(Equation):

    def __init__(self, var, terminationValue = 1e10):
        """

        The `var` argument must contain both positive and negative to
        define the zero level set.

        The `terminationValue` represents the furthest distance from
        the interface that the signed distance function will be
        calculated.
        
        """
        
        self.mesh = var.getMesh()
	self.terminationValue = terminationValue
        
	Equation.__init__(
            self,
            var = var,
            terms = (),
            solver = None)

        self._calcNumericQuantities()

    def solve(self):
        
        setValueFlag = self._calcInterfaceValues()
        setValueFlag = self._calcInitialTrialValues(setValueFlag)
        self._calcRemainingValues(setValueFlag)
        
    def _calcNumericQuantities(self):
        self.cellToCellIDs = Numeric.array(self.mesh.getCellToCellIDsFilled())
        self.cellToCellDistances = Numeric.array(self.mesh.getCellToCellDistances())
        self.cellNormals = Numeric.array(self.mesh.getCellNormals())
        self.cellAreas = Numeric.array(self.mesh.getCellAreas())

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
           >>> v = 0.5 * Numeric.sqrt(2) / 2
           >>> answer = Numeric.array((-1,-0.5,-1,-0.5, v,-0.5,-1,-0.5,-1))
           >>> Numeric.allclose(answer, Numeric.array(var))
           1

           >>> mesh = Grid2D(dx = .5, dy = 2., nx = 3, ny = 3)
           >>> from distanceVariable import DistanceVariable
           >>> var = DistanceVariable(mesh = mesh, value =
           ...     (-1, -1, -1, -1, 1, -1, -1, -1, -1))
           >>> DistanceEquation(var)._calcInterfaceValues()
           [0,1,0,1,1,1,0,1,0,]
           >>> v = 0.25 / Numeric.sqrt(2)
           >>> answer = Numeric.array((-1,-1,-1,-0.25,v,-0.25,-1,-1,-1))
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

        distance = MA.masked_values(abs(phi * dAP / (phi - phiAdj)) * cellZeroFlag, 0)

        distance = MA.sort(distance, axis = 1)

        s = distance[:,0]
        t = distance[:,1]

        sign = -1 + 2 * (Numeric.array(self.var) > 0)

        signedDistance = MA.where(s.mask(),
                                 self.var,
                                 MA.where(t.mask(),
                                          sign * s,
                                          sign * s * t / MA.sqrt(s**2 + t**2)))

##        argmins = MA.argmin(abs(distance), axis = 1)

##        argmins = Numeric.transpose(Numeric.array((Numeric.arange(N),argmins)))

##        argmins = argmins[:,0] * M + argmins[:,1]

##        print 'distance before',distance

##        distance = MA.take(distance.flat, argmins)

##        print 'distance',distance
        
        cellFlag = Numeric.sum(cellZeroFlag, axis = 1)

##        self.var.setValue(MA.where(cellFlag > 0, distance, self.var))
        self.var.setValue(signedDistance)

        return Numeric.where(cellFlag > 0, 1, 0)
    
    def _calcLinear(self, phi, d, cellID, adjCellID):
        return (phi + d, phi - d)
    
    def _calcQuadratic(self, phi1, phi2, n1, n2, d1, d2, area1, area2, cellID, adjCellID1, adjCellID2):
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
           >>> Numeric.allclose(answer, eqn._calcQuadratic(1, 1, n1, n2, 1, 1, 0, 0, 0, 0, 0))
           1
           >>> n1 = Numeric.array((-1,-1)) / sqrt
           >>> n2 = Numeric.array((1,-1)) / sqrt
           >>> d = 1
           >>> answer = Numeric.array((1 + d, 1 - d))
           >>> Numeric.allclose(answer, eqn._calcQuadratic(1, 1, n1, n2, d * sqrt, d * sqrt, 0, 0, 0, 0, 0))
           1
           >>> n1 = Numeric.array((0,-1))
           >>> n2 = Numeric.array((-1,-1)) / sqrt
           >>> d = 0.3
           >>> answer = Numeric.array((1 + d, 1 - d))
           >>> Numeric.allclose(answer, eqn._calcQuadratic(1, 1, n1, n2, d, d * sqrt, 0, 0, 0, 0, 0))
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
           >>> Numeric.allclose(Numeric.array((-.5 / sqrt, .5, .5, .5 + 1 / sqrt)), Numeric.array(var))
           1

           >>> mesh = Grid2D(dx = .5, dy = 2., nx = 3, ny = 3)
           >>> var = DistanceVariable(mesh = mesh, value = (-1, -1, -1, -1, 1, -1, -1, -1, -1))
           >>> eqn = DistanceEquation(var)
           >>> setValueFlag = eqn._calcInterfaceValues()
           >>> eqn._calcInitialTrialValues(setValueFlag)
           [2,1,2,1,1,1,2,1,2,]
           >>> v = -1.40771446
           >>> v1 = 0.25 / Numeric.sqrt(2)
           >>> Numeric.allclose(Numeric.array((v, -1, v, -.25, v1, -.25, v, -1, v)), Numeric.array(var), atol = 1e-6)
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

        trialCellIDs = Numeric.nonzero(Numeric.where(setValueFlag == 2, 1, 0))

        for cellID in trialCellIDs:
            self.var[cellID] = self._calcTrialValue(cellID, setValueFlag)

        return setValueFlag

    def _calcTrialValue(self, cellID, setValueFlag):

        cellToCellIDs = self.cellToCellIDs[cellID]
        dAP = self.cellToCellDistances[cellID]

        phiAdj = Numeric.take(self.var, cellToCellIDs)

        mask = (Numeric.take(setValueFlag, cellToCellIDs)!=1)
        values = Numeric.where(mask, 1e10, phiAdj)
        
        argsort = Numeric.argsort(Numeric.absolute(values))

        values = Numeric.take(values, argsort)
        dAP = Numeric.take(self.cellToCellDistances[cellID], argsort)
        normals = Numeric.take(self.cellNormals[cellID], argsort)
        areas = Numeric.take(self.cellAreas[cellID], argsort)

        sign = -1
        if self.var[cellID] > 0:
            sign = 1

        NSetValues = len(mask) - Numeric.sum(mask)

        
        if NSetValues == 0:
            raise Error
        elif NSetValues == 1:
            lin = self._calcLinear(values[0], dAP[0], cellID, cellToCellIDs[argsort[0]])
            if sign > 0:
                return lin[0]
            else:
                return lin[1]
        elif NSetValues >= 2:
            quad = self._calcQuadratic(values[0], values[1], normals[0], normals[1], dAP[0], dAP[1], areas[0], areas[1], cellID, cellToCellIDs[argsort[0]], cellToCellIDs[argsort[1]])
            if sign > 0:
                return quad[0]
            else:
                return quad[1]

    def _calcRemainingValues(self, setValueFlag):
        """
        
        After the initial trial values have been evaluated this
        routine marches out to evaluate all the remaining cells.
        
        >>> from fipy.meshes.grid2D import Grid2D
        >>> mesh = Grid2D(dx = 1., dy = 1., nx = 1, ny = 4)
        >>> from distanceVariable import DistanceVariable
        >>> var = DistanceVariable(mesh = mesh, value = (-1, 1, 1, 1))
        >>> eqn = DistanceEquation(var)
        >>> setValueFlag = eqn._calcInterfaceValues()
        >>> setValueFlag = eqn._calcInitialTrialValues(setValueFlag)
        >>> eqn._calcRemainingValues(setValueFlag)
        >>> Numeric.allclose(Numeric.array((-.5, .5, 1.5, 2.5)), Numeric.array(var))
        1

        >>> mesh = Grid2D(dx = 1., dy = 1., nx = 3, ny = 3)
        >>> var = DistanceVariable(mesh = mesh, value = (-1, 1, 1, 1, 1, 1, 1, 1, 1))
        >>> eqn = DistanceEquation(var)
        >>> setValueFlag = eqn._calcInterfaceValues()
        >>> setValueFlag = eqn._calcInitialTrialValues(setValueFlag)
        >>> eqn._calcRemainingValues(setValueFlag)
        >>> x = mesh.getCellCenters()[:,0]
        >>> y = mesh.getCellCenters()[:,1]
        >>> answer = Numeric.sqrt((x - .5)**2 + (y - .5)**2) - 0.5
        >>> Numeric.allclose(answer, Numeric.array(var), atol = 5e-1)
        1
        
        """
        trialCellIDs = list(Numeric.nonzero(Numeric.where(setValueFlag == 2, 1, 0)))

        values = Numeric.take(self.var, trialCellIDs)
        argmin = Numeric.argmin(abs(values))
        cellID = trialCellIDs[argmin]

        while len(trialCellIDs) > 0 and abs(self.var[cellID]) < self.terminationValue:

            values = Numeric.take(self.var, trialCellIDs)
            argmin = Numeric.argmin(abs(values))
            cellID = trialCellIDs[argmin]
            setValueFlag[cellID] = 1
            trialCellIDs.remove(cellID)
            cellToCellIDs = self.mesh.getCellToCellIDsFilled()[cellID]
            for adjCellID in cellToCellIDs:
                if setValueFlag[adjCellID] != 1:
                    self.var[adjCellID] = self._calcTrialValue(adjCellID, setValueFlag)
                    if setValueFlag[adjCellID] == 0:
                        setValueFlag[adjCellID] = 2
                        trialCellIDs.append(adjCellID)

                        
            

            
def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__": 
    _test() 




