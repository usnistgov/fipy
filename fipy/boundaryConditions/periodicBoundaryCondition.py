#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - a finite volume PDE solver in Python
 # 
 #  FILE: "periodicBoundaryCondition.py"
 #                                    created: 11/18/04 {4:31:51 PM} 
 #                                last update: 11/20/04 {11:18:18 PM} 
 #  Author: Jonathan Guyer
 #  E-mail: guyer@nist.gov
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #  
 # ========================================================================
 # This document was prepared at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this document is not subject to copyright
 # protection and is in the public domain.  periodicBoundaryCondition.py
 # is an experimental work.  NIST assumes no responsibility whatsoever
 # for its use by other parties, and makes no guarantees, expressed
 # or implied, about its quality, reliability, or any other characteristic.
 # We would appreciate acknowledgement if the document is used.
 # 
 # This document can be redistributed and/or modified freely
 # provided that any derivative works bear some notice that they are
 # derived from it, and any modified versions bear some notice that
 # they have been modified.
 # ========================================================================
 #  See the file "license.terms" for information on usage and  redistribution
 #  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #  
 #  Description: 
 # 
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2004-11-18 JEG 1.0 original
 # ###################################################################
 ##

"""Periodic boundary condition
"""
__docformat__ = 'restructuredtext'

import Numeric

from fipy.boundaryConditions.boundaryCondition import BoundaryCondition
from fipy.tools import array
from fipy.tools.sparseMatrix import SparseMatrix

class PeriodicBoundaryCondition(BoundaryCondition):
    def __init__(self, faces1, faces2):
	self.faces1 = faces1
	self.faces2 = faces2
	
	if len(faces1) is not len(faces2):
	    raise "Incompatible numbers of faces"
	
	self.face1Ids = Numeric.array([face.getID() for face in self.faces1])
	self.face2Ids = Numeric.array([face.getID() for face in self.faces2])
	self.adjacentCell1Ids = Numeric.array([face.getCellID() for face in self.faces1])
	self.adjacentCell2Ids = Numeric.array([face.getCellID() for face in self.faces2])

    def buildMatrix(self, Ncells, MaxFaces, cell1dia, cell1off, coeffScale):
	"""Modify **L** to make `faces1` and `faces2` contiguous.
	**b** is unchanged.
	
	A `tuple` of (`LL`, `bb`) is calculated, to be added to the 
	Term's (**L**, **b**) matrices.
	
	:Parameters:
	    
	  - `Ncells`:   Size of matrices
	  - `MaxFaces`: bandwidth of **L**
	  - `cell1dia`: contribution to adjacent cell diagonal by this 
	    exterior face	    
	  - `cell1off`: contribution to this cell diagonal by adjacent 
	    exterior face
	  - `coeffScale`: dimensionality of the coefficient matrix
	"""
	
	LL = SparseMatrix(size = Ncells, bandwidth = MaxFaces)
	
	diagonalContribution = array.take(cell1dia[:],self.face1Ids) / (2 * coeffScale)
	offdiagonalContribution = array.take(cell1off[:],self.face2Ids) / (2 * coeffScale)
	
	LL.addAt(diagonalContribution, self.adjacentCell1Ids, self.adjacentCell1Ids)
	LL.addAt(offdiagonalContribution, self.adjacentCell1Ids, self.adjacentCell2Ids)
	LL.addAt(offdiagonalContribution, self.adjacentCell2Ids, self.adjacentCell1Ids)
	LL.addAt(diagonalContribution, self.adjacentCell2Ids, self.adjacentCell2Ids)
	
	return (LL, 0)
	
    def getDerivative(self, order):
	return self, self



