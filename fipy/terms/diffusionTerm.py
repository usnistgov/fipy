"""
-*-Pyth-*-
###################################################################
 Alpha - Core code development for Alpha

 FILE: "diffusionTerm.py"
                                   created: 11/13/03 {11:37:00 AM} 
                               last update: 11/14/03 {5:07:34 PM} 
 Author: Jonathan Guyer
 E-mail: jguyer@his.com
   mail: Alpha Cabal
    www: http://alphatcl.sourceforge.net
 
Copyright (c) 2003  Jonathan Guyer

See the file "license.terms" for information on usage and redistribution
of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 See the file "license.terms" for information on usage and  redistribution
 of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 
###################################################################
"""

import faceTerm

class DiffusionTerm(faceTerm.FaceTerm):
    def __init__(self,equation,diffCoeff):
        stencil = (1,1)
	faceTerm.FaceTerm.__init__(self, stencil , equation)
	self.diffCoeff = diffCoeff
	
    def updateCoeff(self,dt):
	faces = self.equation.mesh().faces()
	self.coeff = Numeric.zeroes(len(faces),'d')
	for face in faces:
	    self.coeff[face.id()] = self.diffCoeff * face.area() / face.cellDistance()
	
	


