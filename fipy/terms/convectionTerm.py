"""
-*-Pyth-*-
###################################################################
 PFM - Python-based phase field solver

 FILE: "convectionTerm.py"
                                   created: 11/13/03 {11:39:03 AM} 
                               last update: 12/3/03 {10:37:09 AM} 
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
 
 Description: 

 History

 modified   by  rev reason
 ---------- --- --- -----------
 2003-11-13 JEG 1.0 original
###################################################################
"""

from faceTerm import FaceTerm
import Numeric

class ConvectionTerm(FaceTerm):
    def __init__(self, convCoeff, mesh, boundaryConditions, stencil):
	stencil = ((0.5,-0.5), 'None')
	FaceTerm.__init__(self,stencil,mesh,boundaryConditions)
	self.convCoeff = convCoeff
	
    def updateCoeff(self,dt):
	areas = self.mesh.getOrientedAreaProjections()
	self.coeff = meshes.tools.arraySqrtDot(self.convCoeff,areas)

    def setConvCoeff(self,convCoeff):
        """
        Added to allow the convection coefficient to be updated
        in the equation module
        """
        self.convCoeff = convCoeff
	


