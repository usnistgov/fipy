"""
-*-Pyth-*-
###################################################################
 PFM - Python-based phase field solver

 FILE: "variable.py"
                                   created: 11/10/03 {3:15:38 PM} 
                               last update: 12/9/03 {2:39:10 PM} 
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
 2003-11-10 JEG 1.0 original
###################################################################
"""

import meshes.tools
import Numeric

class Variable:
    
    def __init__(self, mesh, name = '', value=0.):
	self.mesh = mesh
	self.name = name
	self.value = value
	    
    def __coerce__(self,other):
	if type(other) in [type(Numeric.array((1.))),type(1.0),type(1),type(True)]:
	    return (self.getValue(),other)
	else:
	    return None
	
    def __getitem__(self, index): 
	return self.getValue()[index]

    def __setitem__(self, index, value): 
	self.value[index] = value
		
    def getMesh(self):
	return self.mesh

    def getValue(self):
        return self.value

