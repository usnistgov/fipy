## -*-Pyth-*-
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "faceTerm.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 11/21/03 {11:48:08 AM} 
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

class FaceTerm(Term):
    def __init__(self,stencil,faces,interiorFaces,boundaryConditions):
	"""stencil = [phi_adj, phi]
	"""
	Term.__init__(self,stencil)
        self.faces = faces
        self.interiorFaces = interiorFaces
        self.boundaryConditions = boundaryConditions
        
    def buildMatrix(self,L,array,b):
	for face in self.interiorFaces:
            id1 = face.getCellId(0)
            id2 = face.getCellId(1)
            faceId = face.getId()
            L[id1,id1]+=self.coeff[faceId] * self.stencil[1]
            L[id1,id2]-=self.coeff[faceId] * self.stencil[0]
            L[id2,id1]-=self.coeff[faceId] * self.stencil[0]
            L[id2,id2]+=self.coeff[faceId] * self.stencil[1]

        for boundaryCondition in self.boundaryConditions:
            for face in boundaryCondition.getFaces():
                cellId = face.getCellId()
                faceId = face.getId()
                LL,bb = boundaryCondition.update(face,self.coeff[faceId],self.stencil)
                L[cellId,cellId] += LL
                b[cellId] += bb
