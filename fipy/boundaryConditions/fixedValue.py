#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "fixedValue.py"
 #                                    created: 11/15/03 {9:47:59 PM} 
 #                                last update: 12/15/04 {4:52:11 PM}
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #  Author: James Warren <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
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
__docformat__ = 'restructuredtext'

import Numeric

from fipy.boundaryConditions.boundaryCondition import BoundaryCondition
from fipy.tools import array
from fipy.tools import vector
from fipy.tools.sparseMatrix import SparseMatrix

class FixedValue(BoundaryCondition):
    def buildMatrix(self, Ncells, MaxFaces, coeff):
	"""Set boundary equal to value.
	
	A `tuple` of (`LL`, `bb`) is calculated, to be added to the 
	Term's (**L**, **b**) matrices.
	
	:Parameters:
	    
	  - `Ncells`:   Size of matrices
	  - `MaxFaces`: bandwidth of **L**
	  - `coeff`:    contribution to adjacent cell diagonal and **b**-vector by 
	                this exterior face
	"""
	
	LL = SparseMatrix(size = Ncells, bandwidth = MaxFaces)
	LL.addAt(array.take(coeff['cell 1 diag'],self.faceIds), self.adjacentCellIds, self.adjacentCellIds)
	
	bb = Numeric.zeros((Ncells,),'d')
	vector.putAdd(bb, self.adjacentCellIds, array.take(-coeff['cell 1 offdiag'],self.faceIds) * self.getValue())

	return (LL, bb)
	
    def getValue(self):
	return self.value



