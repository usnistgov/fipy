#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  PyFiVol - Python-based finite volume PDE solver
 # 
 #  FILE: "faceTerm.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 2/18/04 {11:37:27 AM} 
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

from fivol.terms.term import Term
import fivol.tools.vector
import fivol.tools.array as array
from fivol.inline import inline

class FaceTerm(Term):
    def __init__(self,weight,mesh,boundaryConditions):
	Term.__init__(self,weight)
        self.mesh = mesh
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
    def implicitBuildMatrix(self, L, coeffScale, id1, id2, b, varScale):
        L.update_add_pyarray_at_indices(self.implicit['cell 1 diag'][:self.interiorN] / coeffScale,id1,id1)
        L.update_add_pyarray_at_indices(self.implicit['cell 1 offdiag'][:self.interiorN] / coeffScale,id1,id2)
        L.update_add_pyarray_at_indices(self.implicit['cell 2 offdiag'][:self.interiorN] / coeffScale,id2,id1)
        L.update_add_pyarray_at_indices(self.implicit['cell 2 diag'][:self.interiorN] / coeffScale,id2,id2)

        for boundaryCondition in self.boundaryConditions:
            LL,bb,ids = boundaryCondition.getContribution(self.implicit['cell 1 diag'],self.implicit['cell 1 offdiag'])
                
            L.update_add_pyarray_at_indices(LL/coeffScale,ids,ids)
            ## WARNING: the next line will not work if one cell has two faces on the same
            ## boundary. Numeric.put will not add both values to the b array but over write
            ## the first with the second. We really need a putAdd function rather than put.
            ## Numeric.put(b,ids,Numeric.take(b,ids)+bb)
		
            fivol.tools.vector.putAdd(b, ids, bb/(coeffScale * varScale))

    def explicitBuildMatrix(self, oldArray, id1, id2, b, coeffScale, varScale):

        inline.optionalInline(self._explicitBuildMatrixIn, self._explicitBuildMatrixPy, oldArray, id1, id2, b, coeffScale, varScale)
        
        for boundaryCondition in self.boundaryConditions:

            LL,bb,ids = boundaryCondition.getContribution(self.explicit['cell 1 diag'],self.explicit['cell 1 offdiag'])
            oldArrayIds = array.take(oldArray, ids)
            fivol.tools.vector.putAdd(b, ids, -LL * oldArrayIds/(coeffScale * varScale))
            fivol.tools.vector.putAdd(b, ids, bb/(coeffScale * varScale))

    def _explicitBuildMatrixIn(self, oldArray, id1, id2, b, coeffScale, varScale):
        try:
            cell1Diag = self.explicit['cell 1 diag'].getNumericValue()
            cell1OffDiag = self.explicit['cell 1 offdiag'].getNumericValue()
            cell2Diag = self.explicit['cell 2 diag'].getNumericValue()
            cell2OffDiag = self.explicit['cell 2 offdiag'].getNumericValue()
        except:
            cell1Diag = self.explicit['cell 1 diag'][:]
            cell1OffDiag = self.explicit['cell 1 offdiag'][:]
            cell2Diag = self.explicit['cell 2 diag'][:]
            cell2OffDiag = self.explicit['cell 2 offdiag'][:]

        weight = self.weight['explicit']
        try:
            coeff = self.coeff.getNumericValue()
        except:
            coeff = self.coeff
            
        inline.runInlineLoop1("""
            double oldArrayId1 = oldArray(id1(i));
            double oldArrayId2 = oldArray(id2(i));

            b(id1(i)) += -coeff(i) * (cell1Diag * oldArrayId1 + cell1OffDiag * oldArrayId2) / coeffScale;
            b(id2(i)) += -coeff(i) * (cell2Diag * oldArrayId2 + cell2OffDiag * oldArrayId1) / coeffScale;
        """,oldArray = oldArray.getNumericValue(),
            id1 = id1,
            id2 = id2,
            b = b,
            coeffScale = coeffScale,
            cell1Diag = weight['cell 1 diag'],
            cell1OffDiag = weight['cell 1 offdiag'],
            cell2Diag = weight['cell 2 diag'],
            cell2OffDiag = weight['cell 2 offdiag'],
            coeff = coeff,
            ni = self.interiorN)
    
    def _explicitBuildMatrixPy(self, oldArray, id1, id2, b, coeffScale, varScale):

        oldArrayId1 = array.take(oldArray, id1)
        oldArrayId2 = array.take(oldArray, id2)

        fivol.tools.vector.putAdd(b, id1, -(self.explicit['cell 1 diag'][:self.interiorN] * oldArrayId1[:] + self.explicit['cell 1 offdiag'][:self.interiorN] * oldArrayId2[:])/coeffScale)
        fivol.tools.vector.putAdd(b, id2, -(self.explicit['cell 2 diag'][:self.interiorN] * oldArrayId2[:] + self.explicit['cell 2 offdiag'][:self.interiorN] * oldArrayId1[:])/coeffScale)

                 
    def buildMatrix(self,L,oldArray,b,coeffScale,varScale):
	"""Implicit portion considers
	"""
	
	id1, id2 = self.mesh.getAdjacentCellIDs()
	id1 = id1[:self.interiorN]
	id2 = id2[:self.interiorN]
	
        ## implicit
        if self.weight.has_key('implicit'):
	    self.implicitBuildMatrix(L, coeffScale, id1, id2, b, varScale)

        if self.weight.has_key('explicit'):
            self.explicitBuildMatrix(oldArray, id1, id2, b, coeffScale, varScale)

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
		
##                fivol.tools.vector.putAdd(b, ids, bb/(coeffScale * varScale))

##        if self.weight.has_key('explicit'):

##	    oldArrayId1 = array.take(oldArray, id1)
##	    oldArrayId2 = array.take(oldArray, id2)

##            fivol.tools.vector.putAdd(b, id1, -(self.explicit['cell 1 diag'][:self.interiorN] * oldArrayId1[:] + self.explicit['cell 1 offdiag'][:self.interiorN] * oldArrayId2[:])/coeffScale)
##            fivol.tools.vector.putAdd(b, id2, -(self.explicit['cell 2 diag'][:self.interiorN] * oldArrayId2[:] + self.explicit['cell 2 offdiag'][:self.interiorN] * oldArrayId1[:])/coeffScale)

##            for boundaryCondition in self.boundaryConditions:

##                LL,bb,ids = boundaryCondition.getContribution(self.explicit['cell 1 diag'],self.explicit['cell 1 offdiag'])
##                oldArrayIds = array.take(oldArray, ids)
##                fivol.tools.vector.putAdd(b, ids, -LL * oldArrayIds/(coeffScale * varScale))
##                fivol.tools.vector.putAdd(b, ids, bb/(coeffScale * varScale))
