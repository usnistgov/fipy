#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "vectorFaceVariable.py"
 #                                    created: 12/9/03 {3:22:07 PM} 
 #                                last update: 9/3/04 {10:41:45 PM} 
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
 # protection and is in the public domain.  FiPy is an experimental
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

import Numeric

from fipy.variables.variable import Variable

class VectorFaceVariable(Variable):
    def __init__(self,mesh,name = '',value=0., unit = None, length = 1):

	length = len(mesh.getFaces())
	array = Numeric.zeros([length,mesh.getDim()],'d')
	
	Variable.__init__(self, mesh = mesh, name = name, value = value, unit = unit, array = array)

    def getVariableClass(self):
	return VectorFaceVariable

##    def getMag(self):
##	if self.mag is None:
##	    from magVectorVariable import MagVectorVariable
##	    self.mag = MagVectorVariable(self)
	    
##	return self.mag

	


	
