## -*-Pyth-*-
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "faceTerm.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 11/26/03 {10:30:04 AM} 
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

from term import Term
import Numeric

class FaceTerm(Term):
    def __init__(self,stencil,mesh,boundaryConditions):
	"""stencil = [phi_adj, phi]
	"""
	Term.__init__(self,stencil)
        self.mesh = mesh
        self.interiorN = len(self.mesh.getInteriorFaces())
        self.boundaryConditions = boundaryConditions
	
	self.calcAdjacentCellIDs()
	
    def calcAdjacentCellIDs(self):
	self.id1 = Numeric.zeros((self.interiorN))
	self.id2 = Numeric.zeros((self.interiorN))
	faces = self.mesh.getFaces()
	for i in range(self.interiorN):
	    self.id1[i] = faces[i].getCellId(0)
	    self.id2[i] = faces[i].getCellId(1)
	
    def actuallyDoSomething(self, L, aa, bb, id1, id2):
	L.update_add_something(aa,id1,id1)
	L.update_add_something(bb,id1,id2)
	L.update_add_something(bb,id2,id1)
	L.update_add_something(aa,id2,id2)
	
    def buildMatrix(self,L,array,b):
	aa =  self.coeff[:self.interiorN]*self.stencil[1]
	bb = -self.coeff[:self.interiorN]*self.stencil[0]
	
	self.actuallyDoSomething(L, aa, bb, self.id1, self.id2)
	
        for boundaryCondition in self.boundaryConditions:
            for face in boundaryCondition.getFaces():
                cellId = face.getCellId()
                faceId = face.getId()
                LL,bb = boundaryCondition.update(face,self.coeff[faceId],self.stencil)
                L[cellId,cellId] += LL
                b[cellId] += bb
