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
        positiveIDs = Set([])

        ## obtain positive and negative IDs

        for id in unevaluatedIDs:
            if self.value[id] > 0:
                positiveIDs.add(id)

        negativeIDs = unevaluatedIDs - positiveIDs

        ## obtain interface IDs

        positiveNeighborIDs = Set([])
        negativeNeighborIDs = Set([])

        for id in positiveIDs:
            positiveNeighborIDs.union_update(Set(self.mesh.getAdjacentCellIDs()[id]))

        for id in negativeIDs:
            negativeNeighborIDs.union_update(Set(self.mesh.getAdjacentCellIDs()[id]))

        interfaceIDs = positiveNeighborIDs & negativeNeighborIDs

        ## calculate interface values

        tmpValue = self.value.copy()

        for id in interfaceIDs:
            tmpValue[id] = self.setInterfaceValue(id)

        self.value = tmpValue

        ## find trial IDs

        evaluatedIDs = interfaceIDs
        unevaluatedIDs = unevaluatedIDs - interfaceIDs
        trialIDs = Set[()]
        
        for id in evaluatedIDs:
            adjIDs = Set(self.mesh.adjacentCellIDs()[id])
            trialCellIDs.union_update(adjIDs & unevaluatedIDs)

        ## calculate trial cell IDs

        for id in trialCellIDs:
            self.value[id] = self.setTrialValue(id)
        
        unevaluatedIDs = unevaluatedIDs - trialIDs

        ## calculate remaining unevaluted IDs
        
        while len(trialIDs) > 0:

            id = self.getMinimumID(trialIDs)

            evaluatedIds.add(id)
            trialIDs.remove(id)

            adjIDs = Set(self.mesh.adjacentCellIDs()[id])
            newTrialIDs = adjIDs & unevaluatedIDs

            for trialID in newTrialIDs:
                self.setTrialValue(trialID)

            trialsIDs.union_update(newTrialsIDs)

            if abs(self.value[id]) > self.narrowBandWidth / 2:
                break
            

    def getMinimumID(self, set):
        IDs = Numeric.array(set)
        vals = Numeric.take(self.value, IDs)
        arg = Numeric.argsort(vals)
        return = IDs[arg]

    def setTrialValue(self):
        pass
        
    def setInterfaceValue(self):
        pass
