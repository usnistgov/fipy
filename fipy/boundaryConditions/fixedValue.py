#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  PyFiVol - Python-based finite volume PDE solver
 # 
 #  FILE: "fixedValue.py"
 #                                    created: 11/15/03 {9:47:59 PM} 
 #                                last update: 4/2/04 {4:05:53 PM}
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

"""Fixed value (Dirichlet) boundary condition
"""

from fipy.boundaryConditions.boundaryCondition import BoundaryCondition
from fipy.tools import array

class FixedValue(BoundaryCondition):
    """Fixed value (Dirichlet) boundary condition
    """
    
    def update(self,face,cell1dia,cell1off):
	"""???
	
	Arguments:
	    
	    'face' -- *unused*
	    
	    'cell1dia' -- contribution to adjacent cell diagonal by this exterior face
	    
	    'cell1off' -- contribution to b-vector by this exterior face
	"""
        return (cell1dia,-cell1off*self.value)

    def getContribution(self,cell1dia,cell1off):
	"""???
	
	Arguments:
	    
	    'cell1dia' -- contribution to adjacent cell diagonal by this exterior face
	    
	    'cell1off' -- contribution to b-vector by this exterior face
	"""
	return (array.take(cell1dia[:],self.faceIds),
		array.take(-cell1off[:],self.faceIds)*self.value,
		self.adjacentCellIds)



