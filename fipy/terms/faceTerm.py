#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "faceTerm.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 7/15/04 {3:34:26 PM} 
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
 # system.  NIST assumes no responsibility whatsever for its use by
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
 #  2003-11-17 JEG 1.0 original
 # ###################################################################
 ##
 
import Numeric

from fipy.terms.term import Term
import fipy.tools.vector
import fipy.tools.array as array
from fipy.tools.inline import inline
from fipy.tools.sparseMatrix import SparseMatrix

from fipy.viewers.gist1DViewer import Gist1DViewer

class FaceTerm(Term):
    def __init__(self,weight,mesh,boundaryConditions):
	Term.__init__(self, mesh = mesh, weight = weight)
        self.interiorN = len(self.mesh.getInteriorFaces())
        self.boundaryConditions = boundaryConditions

	if self.weight.has_key('implicit'):
	    weight = self.weight['implicit']
	    self.implicit = {
		'cell 1 diag': self.coeff * weight['cell 1 diag'],
		'cell 1 offdiag': self.coeff * weight['cell 1 offdiag'],
		'cell 2 diag': self.coeff * weight['cell 2 diag'],
		'cell 2 offdiag': self.coeff * weight['cell 2 offdiag']
	    }

	if self.weight.has_key('explicit'):
	    weight = self.weight['explicit']
	    self.explicit = {
		'cell 1 diag': self.coeff * weight['cell 1 diag'],
		'cell 1 offdiag': self.coeff * weight['cell 1 offdiag'],
		'cell 2 diag': self.coeff * weight['cell 2 diag'],
		'cell 2 offdiag': self.coeff * weight['cell 2 offdiag']
	    }
	    
## 	import fipy.terms.convectionTerm
## 	if isinstance(self, fipy.terms.convectionTerm.ConvectionTerm):
## 	    self.coeffViewer = Gist1DViewer(vars = (self.coeff,), title = "stupid", minVal = -1000, maxVal = 1000)
            
    def implicitBuildMatrix(self, L, coeffScale, id1, id2, b, varScale):

	L.addAt(array.take(self.implicit['cell 1 diag'][:], self.mesh.getInteriorFaceIDs()) / coeffScale,id1,id1)
	L.addAt(array.take(self.implicit['cell 1 offdiag'][:], self.mesh.getInteriorFaceIDs()) / coeffScale,id1,id2)
	L.addAt(array.take(self.implicit['cell 2 offdiag'][:], self.mesh.getInteriorFaceIDs()) / coeffScale,id2,id1)
	L.addAt(array.take(self.implicit['cell 2 diag'][:], self.mesh.getInteriorFaceIDs()) / coeffScale,id2,id2)
	
        for boundaryCondition in self.boundaryConditions:
            LL,bb,ids = boundaryCondition.getContribution(self.implicit['cell 1 diag'],self.implicit['cell 1 offdiag'])
                
	    L.addAt(LL / coeffScale,ids,ids)
            ## WARNING: the next line will not work if one cell has two faces on the same
            ## boundary. Numeric.put will not add both values to the b array but over write
            ## the first with the second. We really need a putAdd function rather than put.
            ## Numeric.put(b,ids,Numeric.take(b,ids)+bb)
		
            fipy.tools.vector.putAdd(b, ids, bb/(coeffScale * varScale))

    def explicitBuildMatrix(self, oldArray, id1, id2, b, coeffScale, varScale):

        inline.optionalInline(self._explicitBuildMatrixIn, self._explicitBuildMatrixPy, oldArray, id1, id2, b, coeffScale, varScale)
        
        for boundaryCondition in self.boundaryConditions:

            LL,bb,ids = boundaryCondition.getContribution(self.explicit['cell 1 diag'],self.explicit['cell 1 offdiag'])
            oldArrayIds = array.take(oldArray, ids)
            fipy.tools.vector.putAdd(b, ids, -LL * oldArrayIds/(coeffScale * varScale))
            fipy.tools.vector.putAdd(b, ids, bb/(coeffScale * varScale))

    def _explicitBuildMatrixIn(self, oldArray, id1, id2, b, coeffScale, varScale):

        weight = self.weight['explicit']
        coeff = Numeric.array(self.coeff)
        Nfac = self.mesh.getNumberOfFaces()

        cell1Diag = Numeric.resize(Numeric.array(weight['cell 1 diag']), (Nfac,))
        cell1OffDiag = Numeric.resize(Numeric.array(weight['cell 1 offdiag']), (Nfac,))
        cell2Diag = Numeric.resize(Numeric.array(weight['cell 2 diag']), (Nfac,))
        cell2OffDiag = Numeric.resize(Numeric.array(weight['cell 2 offdiag']), (Nfac,))

	inline.runInlineLoop1("""
	    long int faceID = faceIDs(i);
	    long int cellID1 = id1(i);
	    long int cellID2 = id2(i);
	    double oldArrayId1 = oldArray(cellID1);
	    double oldArrayId2 = oldArray(cellID2);
	 
	    b(cellID1) += -coeff(faceID) * (cell1Diag(faceID) * oldArrayId1 + cell1OffDiag(faceID) * oldArrayId2);
	    b(cellID2) += -coeff(faceID) * (cell2Diag(faceID) * oldArrayId2 + cell2OffDiag(faceID) * oldArrayId1);
	""",oldArray = Numeric.array(oldArray) / coeffScale,
	    id1 = id1,
	    id2 = id2,
	    b = b,
	    cell1Diag = cell1Diag,
	    cell1OffDiag = cell1OffDiag,
	    cell2Diag = cell2Diag,
	    cell2OffDiag = cell2OffDiag,
	    coeff = coeff,
	    faceIDs = self.mesh.getInteriorFaceIDs(),
	    ni = len(self.mesh.getInteriorFaceIDs()))

    def _explicitBuildMatrixPy(self, oldArray, id1, id2, b, coeffScale, varScale):
        oldArrayId1, oldArrayId2 = self.getOldAdjacentValues(oldArray, id1, id2)

	cell1diag = array.take(self.explicit['cell 1 diag'], self.mesh.getInteriorFaceIDs())
	cell1offdiag = array.take(self.explicit['cell 1 offdiag'], self.mesh.getInteriorFaceIDs())
	cell2diag = array.take(self.explicit['cell 2 diag'], self.mesh.getInteriorFaceIDs())
	cell2offdiag = array.take(self.explicit['cell 2 offdiag'], self.mesh.getInteriorFaceIDs())
	
	fipy.tools.vector.putAdd(b, id1, -(cell1diag * oldArrayId1[:] + cell1offdiag * oldArrayId2[:])/coeffScale)
	fipy.tools.vector.putAdd(b, id2, -(cell2diag * oldArrayId2[:] + cell2offdiag * oldArrayId1[:])/coeffScale)

    def getOldAdjacentValues(self, oldArray, id1, id2):
	return array.take(oldArray, id1), array.take(oldArray, id2)

    def buildMatrix(self, oldArray, coeffScale, varScale, dt):
	"""Implicit portion considers
	"""

## 	print self, "coeff:\n", self.coeff
	
## 	import fipy.terms.convectionTerm
## 	if isinstance(self, fipy.terms.convectionTerm.ConvectionTerm):
## 	    self.coeffViewer.plot()

	self.dt = dt
	
	id1, id2 = self.mesh.getAdjacentCellIDs()
	id1 = array.take(id1, self.mesh.getInteriorFaceIDs())
	id2 = array.take(id2, self.mesh.getInteriorFaceIDs())
	
        N = len(oldArray)
        b = Numeric.zeros((N),'d')
        L = SparseMatrix(size = N)
        
        ## implicit
        if self.weight.has_key('implicit'):
	    self.implicitBuildMatrix(L, coeffScale, id1, id2, b, varScale)

        if self.weight.has_key('explicit'):
            self.explicitBuildMatrix(oldArray, id1, id2, b, coeffScale, varScale)
            
        return (L, b)

##    def buildMatrix(self,L,oldArray,b,coeffScale,varScale):
##	"""Implicit portion considers
##	"""
	
##	id1, id2 = self.mesh.getAdjacentCellIDs()
##	id1 = id1[:self.interiorN]
##	id2 = id2[:self.interiorN]
	
##        ## implicit
##        if self.weight.has_key('implicit'):
	    
##            L.update_add_pyarray_at_indices(self.implicit['cell 1 diag'][:self.interiorN] / coeffScale,id1,id1)
##            L.update_add_pyarray_at_indices(self.implicit['cell 1 offdiag'][:self.interiorN] / coeffScale,id1,id2)
##            L.update_add_pyarray_at_indices(self.implicit['cell 2 offdiag'][:self.interiorN] / coeffScale,id2,id1)
##            L.update_add_pyarray_at_indices(self.implicit['cell 2 diag'][:self.interiorN] / coeffScale,id2,id2)
                
##	    for boundaryCondition in self.boundaryConditions:
##		LL,bb,ids = boundaryCondition.getContribution(self.implicit['cell 1 diag'],self.implicit['cell 1 offdiag'])
                
##		L.update_add_pyarray_at_indices(LL/coeffScale,ids,ids)
##                ## WARNING: the next line will not work if one cell has two faces on the same
##                ## boundary. Numeric.put will not add both values to the b array but over write
##                ## the first with the second. We really need a putAdd function rather than put.
##		## Numeric.put(b,ids,Numeric.take(b,ids)+bb)
		
##                fipy.tools.vector.putAdd(b, ids, bb/(coeffScale * varScale))

##        if self.weight.has_key('explicit'):

##	    oldArrayId1 = array.take(oldArray, id1)
##	    oldArrayId2 = array.take(oldArray, id2)

##            fipy.tools.vector.putAdd(b, id1, -(self.explicit['cell 1 diag'][:self.interiorN] * oldArrayId1[:] + self.explicit['cell 1 offdiag'][:self.interiorN] * oldArrayId2[:])/coeffScale)
##            fipy.tools.vector.putAdd(b, id2, -(self.explicit['cell 2 diag'][:self.interiorN] * oldArrayId2[:] + self.explicit['cell 2 offdiag'][:self.interiorN] * oldArrayId1[:])/coeffScale)

##            for boundaryCondition in self.boundaryConditions:

##                LL,bb,ids = boundaryCondition.getContribution(self.explicit['cell 1 diag'],self.explicit['cell 1 offdiag'])
##                oldArrayIds = array.take(oldArray, ids)
##                fipy.tools.vector.putAdd(b, ids, -LL * oldArrayIds/(coeffScale * varScale))
##                fipy.tools.vector.putAdd(b, ids, bb/(coeffScale * varScale))
