#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - a finite volume PDE solver in Python
 # 
 #  FILE: "periodicBoundaryCondition.py"
 #                                    created: 11/18/04 {4:31:51 PM} 
 #                                last update: 3/10/05 {4:38:27 PM} 
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

__docformat__ = 'restructuredtext'

import Numeric

from fipy.boundaryConditions.boundaryCondition import _BoundaryCondition
from fipy.tools import array
from fipy.tools.sparseMatrix import _SparseMatrix

class PeriodicBoundaryCondition(_BoundaryCondition):
    r"""
    
    A `PeriodicBoundaryCondition` creates a periodic boundary across given
    sets of faces.

    .. warning::

        This boundary condition only works for diffusion with a
        gridded mesh. It currently does not work for convection terms
        or the `NthOrderDiffusionTerm`.

    Usage ::

        PeriodicBoundaryCondition(faces1, faces2)

    """
    def __init__(self, faces1, faces2):
        """
        Creates a `PeriodicBoundaryCondition`.

        :Parameters:
          - `faces1` : A `list` or `tuple` of faces.
          - `faces2` : A `list` or `tuple` of faces.

        The faces in each list are matched in order.
        """
        
	if len(faces1) != len(faces2):
	    raise "Incompatible numbers of faces: %d vs %d"%(len(faces1), len(faces2))

	self.faces1 = faces1
	self.faces2 = faces2
	
	self.faceIds = Numeric.array([face.getID() for face in faces1 + faces2])
	self.face1Ids = Numeric.array([face.getID() for face in faces1])
	self.face2Ids = Numeric.array([face.getID() for face in faces2])

	self.adjacentCellIds = Numeric.array([face.getCellID() for face in faces1 + faces2])
	self.offdiagonalCellIds = Numeric.array([face.getCellID() for face in faces2 + faces1])
	
	self.adjacent1CellIds = Numeric.array([face.getCellID() for face in faces1])
	self.adjacent2CellIds = Numeric.array([face.getCellID() for face in faces2])
	
## 	_BoundaryCondition.__init__(self, faces1 + faces2, 0)
	
## 	self.offdiagonalCellIds = Numeric.array([face.getCellID() for face in faces2 + faces1])

    def _buildMatrix(self, Ncells, MaxFaces, coeff):
	"""Modify **L** to make `faces1` and `faces2` contiguous.
	**b** is unchanged.
	
	A `tuple` of (`LL`, `bb`) is calculated, to be added to the 
	Term's (**L**, **b**) matrices.
	
	:Parameters:
	  - `Ncells`:   Size of matrices
	  - `MaxFaces`: bandwidth of **L**
	  - `coeff`: contributions to **L** by this exterior face
	"""
	
	LL = _SparseMatrix(size = Ncells, bandwidth = MaxFaces)

	LL.addAt(array.take(coeff['cell 1 diag'][:] / 2., self.faceIds),self.adjacentCellIds,self.adjacentCellIds)
	LL.addAt(array.take(coeff['cell 1 offdiag'][:] / 2., self.faceIds),self.adjacentCellIds,self.offdiagonalCellIds)
##	LL.addAt(array.take(coeff['cell 2 offdiag'][:], self.faceIds),self.offdiagonalCellIds,self.adjacentCellIds)
##	LL.addAt(array.take(coeff['cell 2 diag'][:], self.faceIds),self.offdiagonalCellIds,self.offdiagonalCellIds)

## 	diagonalContribution = array.take(cell1dia[:],self.faceIds) / (2 * coeffScale)
## 	offdiagonalContribution = array.take(cell1off[:],self.faceIds) / (2 * coeffScale)

## 	LL.addAt(diagonalContribution, self.adjacentCellIds, self.adjacentCellIds)
## 	LL.addAt(offdiagonalContribution, self.adjacentCellIds, self.offdiagonalCellIds)

	return (LL, 0)
	
    def _getDerivative(self, order):
	return self



