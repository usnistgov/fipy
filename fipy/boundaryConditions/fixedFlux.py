#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "fixedFlux.py"
 #                                    created: 11/15/03 {9:47:59 PM} 
 #                                last update: 1/3/07 {2:58:39 PM} 
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

__docformat__ = 'restructuredtext'

from fipy.tools import numerix

from fipy.boundaryConditions.boundaryCondition import BoundaryCondition
from fipy.boundaryConditions.fixedValue import FixedValue
from fipy.tools import vector

class FixedFlux(BoundaryCondition):
    r"""

    The `FixedFlux` boundary condition adds a contribution, equivalent to a
    fixed flux (Neumann condition), to the equation's RHS vector.  The
    contribution, given by `value`, is only added to entries corresponding to
    the specified `faces`, and is weighted by the face areas.
       
    """
    
    def __init__(self,faces,value):
        """
        Creates a `FixedFlux` object.
        
        :Parameters:
            - `faces`: A `list` or `tuple` of `Face` objects to which this condition applies.
            - `value`: The value to impose.
            
        """
        BoundaryCondition.__init__(self,faces,value)
## 	N = len(self.faces)
        ##self.contribution = numerix.zeros((N,),'d')
        # get units right
        ##self.contribution = self.contribution * self.value * self.faces[0].getArea()
        ##for i in range(N):  
        ##    self.contribution[i] = self.value * self.faces[i].getArea()
        self.contribution = self.value * self.faces.getAreas()
        
    def _buildMatrix(self, SparseMatrix, Ncells, MaxFaces, coeff):
        """Leave **L** unchanged and add gradient to **b**
        
        :Parameters:
          - `SparseMatrix`: *unused* (empty matrix)
          - `Ncells`:       Size of **b**-vector
          - `MaxFaces`:     *unused*
          - `coeff`:        *unused*
        """

        bb = numerix.zeros((Ncells,),'d')
## 	vector.putAdd(bb, self.adjacentCellIDs, -numerix.array(self.contribution))
        vector.putAdd(bb, self.adjacentCellIDs, -self.contribution)
        
        return (0, bb)

    def _getDerivative(self, order):
        if order == 1:
            return FixedValue(self.faces, self.value)
        else:
            return BoundaryCondition._getDerivative(self, order)


