#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "addOverFacesVariable.py"
 #                                    created: 4/30/04 {10:39:23 AM} 
 #                                last update: 6/10/04 {4:51:36 PM}
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
 #  2004- 4-30 JEG 1.0 original
 # ###################################################################
 ##

import Numeric

import fipy.tools.array as array
import fipy.tools.inline.inline as inline
from fipy.variables.cellVariable import CellVariable

class AddOverFacesVariable(CellVariable):
    def __init__(self, faceVariable, mesh = None):
        if not mesh:
            mesh = faceVariable.getMesh()

        CellVariable.__init__(self, mesh, hasOld = 0)
    
        self.faceVariable = self.requires(faceVariable)

    def _calcValuePy(self):
        ids = self.mesh.getCellFaceIDs()
        
        contributions = array.take(self.faceVariable[:], ids.flat)

        NCells = self.mesh.getNumberOfCells()

        contributions = array.reshape(contributions,(NCells,-1))
        
        orientations = array.reshape(self.mesh.getCellFaceOrientations(),(NCells,-1))

##         orientations = Numeric.array(orientations)
        
        self.value = array.sum(contributions * orientations,1) / self.mesh.getCellVolumes()
	
    def _calcValueIn(self):

        NCells = self.mesh.getNumberOfCells()
	ids = self.mesh.getCellFaceIDs()

        inline.runInline("""
        int i;
        
        for(i = 0; i < numberOfCells; i++)
          {
          int j;
          value(i) = 0.;
          for(j = 0; j < numberOfCellFaces; j++)
            {
              value(i) += orientations(i,j) * faceVariable(ids(i,j));
            }
            value(i) = value(i) / cellVolume(i);
          }
          """,
              numberOfCellFaces = self.mesh.getMaxFacesPerCell(),
              numberOfCells = NCells,
              faceVariable = self.faceVariable.getNumericValue(),
              ids = Numeric.array(ids),
              value = self.value.value,
              orientations = Numeric.array(self.mesh.getCellFaceOrientations()),
              cellVolume = Numeric.array(self.mesh.getCellVolumes()))
	      

    def calcValue(self):

        inline.optionalInline(self._calcValueIn, self._calcValuePy)



    


