
"""
## -*-Pyth-*-
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "solventVariable.py"
 #                                    created: 12/9/03 {3:02:52 PM} 
 #                                last update: 12/18/03 {12:24:15 AM} 
 #  Author: Jonathan Guyer
 #  E-mail: guyer@nist.gov
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
"""

from componentVariable import ComponentVariable

class SolventVariable(ComponentVariable):
    def __init__(self,mesh,substitutionals):
	ComponentVariable.__init__(self, mesh = mesh, name = "solvent", standardPotential = 0., barrierHeight = 0.)

	self.substitutionals = substitutionals
	for component in self.substitutionals:
	    self.requires(component)

    def calcValue(self):
	self.value[:] = 1.
	for component in self.substitutionals:
	    self.value = self.value - component#.getOld()
