"""
-*-Pyth-*-
###################################################################
 Alpha - Core code development for Alpha

 FILE: "fixedValue.py"
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

import boundaryCondition

class fixedValue(boundaryCondition.BoundaryCondition):
    def __init__(self,faces):
        boundaryCondition.BoundaryCondition.__init__(self,faces)

    def update(self,term):
        for face in self.faces:
            id=face.cells()[0].id()
            term.equation.L()[id1,id1]+=coeff[face.id()] * term.stencil[1]
            term.equation.b()[id1]+=coeff[face.id()] * term.stencil[0] * value
        
    
	


