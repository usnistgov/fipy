#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "solventVariable.py"
 #                                    created: 12/23/03 {4:51:16 PM} 
 #                                last update: 1/15/04 {4:11:39 PM} 
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
 #  See the file "license.terms" for information on usage and  redistribution
 #  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #  
 # ###################################################################
 ##

import Numeric

from componentVariable import ComponentVariable

class SolventVariable(ComponentVariable):
    def __init__(self, mesh, parameters, substitutionals):
	ComponentVariable.__init__(
	    self, 
	    mesh = mesh,
	    parameters = parameters
	    )
	    
	self.concentration = Numeric.array([1])
	for component in substitutionals:
	    self.concentration = self.concentration - component#.getOld()

	self.requires(self.concentration)
	
	self.value = self.concentration[:]
	    
    def calcValue(self):
	self.value = self.concentration[:]
