"""
## -*-Pyth-*-
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "faceTerm.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 12/10/03 {10:11:50 AM} 
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
 
from term import Term
import Numeric

class FaceTerm(Term):
    def __init__(self,weight,mesh,boundaryConditions):
	Term.__init__(self,weight)
        self.mesh = mesh
        self.interiorN = len(self.mesh.getInteriorFaces())
        self.boundaryConditions = boundaryConditions
	
    def buildMatrix(self,L,array,b):
	"""Implicit portion considers
	"""
	
	id1, id2 = self.mesh.getAdjacentCellIDs()
	id1 = id1[:self.interiorN]
	id2 = id2[:self.interiorN]
	
        ## implicit
        if self.weight.has_key('implicit'):
	    weight = self.weight['implicit']
 	    cell1dia = self.coeff*weight['cell 1 diag']
	    cell1off = self.coeff*weight['cell 1 offdiag']
	    cell2dia = self.coeff*weight['cell 2 diag']
	    cell2off = self.coeff*weight['cell 2 offdiag']
	    
	    L.update_add_something(cell1dia[:self.interiorN],id1,id1)
	    L.update_add_something(cell1off[:self.interiorN],id1,id2)
	    L.update_add_something(cell2off[:self.interiorN],id2,id1)
	    L.update_add_something(cell2dia[:self.interiorN],id2,id2)
	    
	    for boundaryCondition in self.boundaryConditions:
		LL,bb,ids = boundaryCondition.getContribution(cell1dia,cell1off)
		L.update_add_something(LL,ids,ids)
		Numeric.put(b,ids,Numeric.take(b,ids)+bb)
		
        ## explicit
        if self.weight.has_key('explicit'):
	    weight = self.weight['explicit']
	    cell1dia = self.coeff*weight['cell 1 diag']
	    cell1off = self.coeff*weight['cell 1 offdiag']
	    cell2off = self.coeff*weight['cell 2 offdiag']
	    cell2dia = self.coeff*weight['cell 2 diag']

            for i in range(self.interiorN):
                b[id1[i]] -= cell1dia[i] * array[id1[i]] + cell1off[i] * array[id2[i]]
                b[id2[i]] -= cell2dia[i] * array[id2[i]] + cell2off[i] * array[id1[i]]
	
            for boundaryCondition in self.boundaryConditions:
                for face in boundaryCondition.getFaces():
                    cellId = face.getCellId()
                    faceId = face.getId()
                    LL,bb = boundaryCondition.update(face,cell1dia[faceId],cell1off[faceId])
                    b[cellId] -= LL * array[cellId]
                    b[cellId] += bb
        
