"""-*-Pyth-*-
###################################################################
 Alpha - Core code development for Alpha

 FILE: "diffusionTerm.py"
                                   created: 11/13/03 {11:37:00 AM} 
                               last update: 11/13/03 {11:43:38 AM} 
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
    def __init__(self,equation,coeff):
	faceTerm.FaceTerm.__init__(self, stencil = (1,1), equation)
	self.coeff = coeff
