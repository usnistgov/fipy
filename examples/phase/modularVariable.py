"""
## -*-Pyth-*-
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "modularVariable.py"
 #                                    created: 12/8/03 {5:47:27 PM} 
 #                                last update: 12/9/03 {2:25:06 PM} 
 #  Author: Jonathan Guyer
 #  E-mail: guyer@nist.gov
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
 #  See the file "license.terms" for information on usage and  redistribution
 #  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #  
 # ###################################################################
 ##
"""

from variables.cellVariable import CellVariable
import Numeric

class ModularVariable(CellVariable):

    def getGrad(self):
        gridSpacing = self.mesh.getMeshSpacing()
        return self.mod(CellVariable.getGrad(self) * gridSpacing) / gridSpacing       

    def getFaceValue(self):
        alpha = self.mesh.getFaceToCellDistanceRatio()
	id1, id2 = self.mesh.getAdjacentCellIDs()
	cell1 = Numeric.take(self[:], id1)
	cell2 = Numeric.take(self[:], id2)
	return self.mod(cell1 - cell2) * alpha + cell2

    def getFaceGrad(self):
        dAP = self.mesh.getCellDistances()
	id1, id2 = self.mesh.getAdjacentCellIDs()
	N = self.mod(Numeric.take(self[:], id2) - Numeric.take(self[:], id1))/dAP
        
	normals = self.mesh.getFaceNormals().copy()
	normals *= Numeric.reshape(self.mesh.getFaceOrientations(),(len(normals),1))
	tangents1 = self.mesh.getFaceTangents1()
	tangents2 = self.mesh.getFaceTangents2()
	cellGrad = self.getGrad()
	grad1 = Numeric.take(cellGrad, id1)
	grad2 = Numeric.take(cellGrad, id2)
	t1grad1 = Numeric.sum(tangents1*grad1,1)
	t1grad2 = Numeric.sum(tangents1*grad2,1)
	t2grad1 = Numeric.sum(tangents2*grad1,1)
	t2grad2 = Numeric.sum(tangents2*grad2,1)
	T1 = (t1grad1 + t1grad2) / 2.
	T2 = (t2grad1 + t2grad2) / 2.
	
	N = Numeric.reshape(N, (len(normals),1)) 
	T1 = Numeric.reshape(T1, (len(normals),1)) 
	T2 = Numeric.reshape(T2, (len(normals),1)) 

	return normals * N + tangents1 * T1 + tangents2 * T2
    
    def mod(self, array):
        pi=Numeric.pi
        return (array + 3. * pi) % (2 * pi) - pi

    def updateOld(self):
        if self.old != None:
            self.value = self.mod(self.value)
	    self.old.setValue(self.value)

    
	


