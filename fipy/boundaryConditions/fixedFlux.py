#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "fixedFlux.py"
 #                                    created: 11/15/03 {9:47:59 PM} 
 #                                last update: 12/8/03 {1:47:41 PM} 
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
 #  2003-11-15 JEG 1.0 original
 # ###################################################################
 ##

"""Fixed flux (Neumann) boundary condition
"""

from boundaryCondition import BoundaryCondition
import Numeric

class FixedFlux(BoundaryCondition):
    """Fixed flux (Neumann) boundary condition
    """
    def __init__(self,faces,value):
	BoundaryCondition.__init__(self,faces,value)
	N = len(self.faces)
	self.contribution = Numeric.zeros((N,),'d')
	for i in range(N):
	    self.contribution[i] = self.value * self.faces[i].getArea()
	
    def update(self,face,cell1dia,cell1off):
	"""Leave L unchanged and add gradient to b
	
	Arguments:
	    
	    'face' -- which 'Face' to update
	    
	    'cell1dia' -- *unused*

	    'cell1off' -- *unused*
	"""
        return (0., self.value * face.getArea())

    def getContribution(self,cell1dia,cell1off):
	"""Leave L unchanged and add gradient to b
	
	Arguments:
	    
	    'cell1dia' -- *unused*

	    'cell1off' -- *unused*
	"""
	return (Numeric.zeros((len(self.faces),),'d'), self.contribution, self.adjacentCellIds)

