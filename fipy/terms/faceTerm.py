"""
-*-Pyth-*-
###################################################################
 PFM - Python-based phase field solver

 FILE: "faceTerm.py"
                                   created: 11/17/03 {10:29:10 AM} 
                               last update: 11/17/03 {3:18:19 PM} 
 Author: Jonathan Guyer
 E-mail: guyer@nist.gov
 Author: Daniel Wheeler
 E-mail: daniel.wheeler@nist.gov
   mail: NIST
    www: http://ctcms.nist.gov
 
========================================================================
This software was developed at the National Institute of Standards
and Technology by employees of the Federal Government in the course
of their official duties.  Pursuant to title 17 Section 105 of the
United States Code this software is not subject to copyright
protection and is in the public domain.  PFM is an experimental
system.  NIST assumes no responsibility whatsoever for its use by
other parties, and makes no guarantees, expressed or implied, about
its quality, reliability, or any other characteristic.  We would
appreciate acknowledgement if the software is used.

This software can be redistributed and/or modified freely
provided that any derivative works bear some notice that they are
derived from it, and any modified versions bear some notice that
they have been modified.
========================================================================
 See the file "license.terms" for information on usage and  redistribution
 of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 
 Description: 

 History

 modified   by  rev reason
 ---------- --- --- -----------
 2003-11-17 JEG 1.0 original
###################################################################
"""

import term

class FaceTerm(term.Term):
    def __init__(self,stencil,equation):
	"""
	stencil = [phi_adj, phi]
	"""
	term.Term.__init__(self,stencil,equation)
	
    def buildMatrix(self):
	var = self.equation.getVar()
	N = len(var)

        mesh = self.equation.getMesh()
	
	for face in mesh.getInteriorFaces():
            cells = face.getCells()
            print face.getId()
            print cells
            id1 = cells[0].getId()
            id2 = cells[1].getId()
            faceId = face.getId()
            self.equation.getL()[id1,id1]+=self.coeff[faceId] * self.stencil[1]
            self.equation.getL()[id1,id2]-=self.coeff[faceId] * self.stencil[0]
            self.equation.getL()[id2,id1]-=self.coeff[faceId] * self.stencil[0]
            self.equation.getL()[id2,id2]+=self.coeff[faceId] * self.stencil[1]

        for boundaryCondition in self.equation.getBoundaryConditions():            
            boundaryCondition.update(self)
