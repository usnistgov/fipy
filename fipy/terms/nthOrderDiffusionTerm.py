#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - a finite volume PDE solver in Python
 # 
 #  FILE: "nthOrderDiffusionTerm.py"
 #                                    created: 5/10/04 {11:24:01 AM} 
 #                                last update: 6/15/04 {11:52:53 AM} 
 #  Author: Jonathan Guyer
 #  E-mail: guyer@nist.gov
 #  Author: Daniel Wheeler
 #  E-mail: daniel.wheeler@nist.gov
 #    mail: NIST
 #     www: http://ctcms.nist.gov
 #  
 # ========================================================================
 # This document was prepared at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this document is not subject to copyright
 # protection and is in the public domain.  FiPy
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
 #  2004-05-10 JEG 1.0 original
 # ###################################################################
 ##

"""

##This `Term` implements a higher order Laplacian term, such as for Cahn-Hilliard

"""
__docformat__ = 'restructuredtext'

import Numeric

from fipy.terms.term import Term
from fipy.tools.sparseMatrix import SparseMatrix
from fipy.tools.sparseMatrix import SparseIdentityMatrix
from fipy.variables.addOverFacesVariable import AddOverFacesVariable
import fipy.tools.array as array
import fipy.tools.vector

class NthOrderDiffusionTerm(Term):
    def __init__(self,coeffs,mesh,boundaryConditions):
	"""
        This term represents a higher order diffusion term. The order of the term is determined
        by the number of `coeffs`, such that
        
	    >>> NthOrderDiffusionTerm(D1, mesh, bcs)
	    
	represents a typical 2nd-order diffusion term of the form
        
        .. raw:: latex
        
           $$ \\nabla\\cdot\\left(D_1 \\nabla \\phi\\right) $$
              
        and
           
            >>> NthOrderDiffusionTerm((D1,D2), mesh, bcs)
            
        represents a 4th-order Cahn-Hilliard term of the form
        
        .. raw:: latex
        
           $$ \\nabla\\cdot\\left[ D_1 \\nabla\\cdot\\left(D_2 \\nabla \\phi\\right) \\right] $$

        and so on.
	"""
	self.order = len(coeffs) * 2
	if len(coeffs) > 0:
	    self.coeff = coeffs[0] * mesh.getFaceAreas() / mesh.getCellDistances()
	else:
	    self.coeff = None
	
	self.boundaryConditions = []
        lowerBoundaryConditions = []
        for bc in boundaryConditions:
            bcDeriv = bc.getDerivative(self.order-2)
	    if bcDeriv:
		self.boundaryConditions.append(bcDeriv)
	    else:
		lowerBoundaryConditions.append(bc)
            
	Term.__init__(self, weight = None, mesh = mesh)
	
	N = mesh.getNumberOfCells()
	if self.order > 0:
	    self.lowerOrderDiffusionTerm = NthOrderDiffusionTerm(coeffs = coeffs[1:], mesh = mesh, boundaryConditions = lowerBoundaryConditions)

    def getCoefficientMatrix(self):
	mesh = self.getMesh()
	
	coefficientMatrix = SparseMatrix(size = mesh.getNumberOfCells(), bandwidth = mesh.getMaxFacesPerCell())
	
	interiorCoeff1 = Numeric.array(self.coeff)
	array.put(interiorCoeff1, mesh.getExteriorFaceIDs(), 0)
        from fipy.variables.faceVariable import FaceVariable
        interiorCoeff1 = FaceVariable(mesh = mesh, value = interiorCoeff1)
        
	contributions = array.take(interiorCoeff1, mesh.getCellFaceIDs())
	contributions = array.sum(contributions, 1)
	
	coefficientMatrix.addAtDiagonal(contributions)
        
	interiorCoeff2 = -array.take(Numeric.array(self.coeff), mesh.getInteriorFaceIDs())
	interiorFaceCellIDs = array.take(mesh.getFaceCellIDs(), mesh.getInteriorFaceIDs())
	
	coefficientMatrix.addAt(interiorCoeff2, interiorFaceCellIDs[:,0], interiorFaceCellIDs[:,1])
	coefficientMatrix.addAt(interiorCoeff2, interiorFaceCellIDs[:,1], interiorFaceCellIDs[:,0])
	
	return coefficientMatrix
	
    def buildMatrix(self, oldArray, coeffScale, varScale):
        if self.order > 0:
            L, b = self.lowerOrderDiffusionTerm.buildMatrix(oldArray, coeffScale, varScale)
            
            coefficientMatrix = self.getCoefficientMatrix()
            self.L = coefficientMatrix * L
            self.b = coefficientMatrix * b
        else:
            N = self.getMesh().getNumberOfCells()
            self.L = SparseIdentityMatrix(size = N)
            self.b = Numeric.zeros((N),'d')

        for boundaryCondition in self.boundaryConditions:
            LL,bb,ids = boundaryCondition.getContribution(self.coeff,-self.coeff)
	    
            self.L.addAt(LL / coeffScale,ids,ids)
            ## WARNING: the next line will not work if one cell has two faces on the same
            ## boundary. Numeric.put will not add both values to the b array but over write
            ## the first with the second. We really need a putAdd function rather than put.
            ## Numeric.put(b,ids,Numeric.take(b,ids)+bb)
                
            fipy.tools.vector.putAdd(self.b, ids, bb/(coeffScale * varScale))
	
	print "order:", self.order
	print "L:", self.L
	print "b:", self.b
	
        return (self.L, self.b)

def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__": 
    _test() 
