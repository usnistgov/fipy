#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "distanceVariable.py"
 #                                    created: 7/29/04 {10:39:23 AM} 
 #                                last update: 10/19/04 {4:38:29 PM} 
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #  
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this software is not subject to copyright
 # protection and is in the public domain.  FiPy is an experimental
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

__docformat__ = 'restructuredtext'

from fipy.variables.cellVariable import CellVariable

class DistanceVariable(CellVariable):
    def __init__(self, mesh = None, name = 'level set variable', value):

        CellVariable(self, mesh, name = name, value = 1.)

    def _calcValue(self):

        Ncells = self.mesh.getNumberOfCells()
        unevaluatedIDs = Set(range(NCells))
        adjCellIDS = self.mesh.getAdjacentCellIDs()
        
        ## obtain positive and negative IDs

        positiveIDs = Set(Numeric.nonzero(self.value > 0))
        negativeIDs = unevaluatedIDs - positiveIDs

        ## obtain interface IDs

        positiveNeighbourIDs = Set(Numeric.take(adjCellIDs, positiveIDs).flat)
        negativeNeighbourIDs = Set(Numeric.take(adjCellIDs, negativeIDs).flat)

        interfaceIDs = positiveNeighborIDs & negativeNeighborIDs
        positiveInterfaceIDs = interfaceIDs & positiveIDs
        negativeInterfaceIDs = interfaceIDs & negativeIDs
        
        ## calculate interface values

        tmpValue = self.value.copy()
        
        for id in positiveInterfaceIDs:
            tmpValue[id] = self.calcInterfaceValue(id, Set(adjCellIDs[id]) & negativeInterfaceIDs, adjCellIDs[id])

        for id in negativeInterfaceIDs:
            tmpValue[id] = self.calcInterfaceValue(id, Set(adjCellIDs[id]) & positiveInterfaceIDs, adjCellIDs[id])

        self.value = tmpValue
        

        ## find trial IDs

        evaluatedIDs = interfaceIDs
        unevaluatedIDs = unevaluatedIDs - evaluatedIDs
        interfaceNeighborIDs = Set(Numeric.take(adjCellIDs, interfaceIDs).flat)
        trialIDs = interfaceNeighbourIDs - interfaceIDs

        ## calculate trial cell IDs

        for id in trialCellIDs:
            self.value[id] = self.calcTrialValue(id, Set(adjCellIDs[id]) & evaluatedIDs)

        unevaluatedIDs = unevaluatedIDs - trialIDs

        ## calculate remaining unevaluted IDs

        trialIDs = list(trialIDs)
         
        while len(trialIDs) > 0:

            id = Numeric.argmin(abs(Numeric.take(self.value, trialsIDs)))
            
            evaluatedIds.add(id)
            trialIDs.remove(id)

            newTrialIDs = Set(adjCellIDs[id]) & unevaluatedIDs

            for trialID in newTrialIDs:
                self.setTrialValue(trialID)
                trialIDs.add(trialID)

            if abs(self.value[id]) > self.narrowBandWidth / 2:
                break
            
    def calcTrialValue(self):
        
        pass
        
    def calcInterfaceValue(self, id, adjIDs, allAdjIDs):
        distances = ()
        val = self.value[adjID]

        for adjID in adjIDs:            
            index = list(allAdjIDs).index(adjID)
            dAP = self.mesh.getCellDistances()[id][index]
            s += (val * dAP / abs(val - self.value[adjID]),)

        for i in len(adjIDs):
            for j in len(adjIDs - 1):
                
