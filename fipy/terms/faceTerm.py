"""
## -*-Pyth-*-
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "faceTerm.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update:  01/05/04 { 2:28:39 PM}
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
 #  2003-11-17 JEG 1.0 original
 # ###################################################################
 ##
"""
 
from term import Term
import Numeric
import tools.vector

class FaceTerm(Term):
    def __init__(self,weight,mesh,boundaryConditions):
	Term.__init__(self,weight)
        self.mesh = mesh
        self.interiorN = len(self.mesh.getInteriorFaces())
        self.boundaryConditions = boundaryConditions
	
	if self.weight.has_key('implicit'):
	    weight = self.weight['implicit']
	    self.implicit = {
		'cell 1 diag': self.coeff*weight['cell 1 diag'],
		'cell 1 offdiag': self.coeff*weight['cell 1 offdiag'],
		'cell 2 diag': self.coeff*weight['cell 2 diag'],
		'cell 2 offdiag': self.coeff*weight['cell 2 offdiag']
	    }
	    
	if self.weight.has_key('explicit'):
	    weight = self.weight['explicit']
	    self.explicit = {
		'cell 1 diag': self.coeff*weight['cell 1 diag'],
		'cell 1 offdiag': self.coeff*weight['cell 1 offdiag'],
		'cell 2 diag': self.coeff*weight['cell 2 diag'],
		'cell 2 offdiag': self.coeff*weight['cell 2 offdiag']
	    }
	    
    def buildMatrix(self,L,oldArray,b):
	"""Implicit portion considers
	"""
	
	id1, id2 = self.mesh.getAdjacentCellIDs()
	id1 = id1[:self.interiorN]
	id2 = id2[:self.interiorN]
	
        ## implicit
        if self.weight.has_key('implicit'):
	    
	    L.update_add_something(self.implicit['cell 1 diag'][:self.interiorN],id1,id1)
	    L.update_add_something(self.implicit['cell 1 offdiag'][:self.interiorN],id1,id2)
	    L.update_add_something(self.implicit['cell 2 offdiag'][:self.interiorN],id2,id1)
	    L.update_add_something(self.implicit['cell 2 diag'][:self.interiorN],id2,id2)
	    
	    for boundaryCondition in self.boundaryConditions:
		LL,bb,ids = boundaryCondition.getContribution(self.implicit['cell 1 diag'],self.implicit['cell 1 offdiag'])
                
		L.update_add_something(LL,ids,ids)
                ## WARNING: the next line will not work if one cell has two faces on the same
                ## boundary. Numeric.put will not add both values to the b array but over write
                ## the first with the second. We really need a putAdd function rather than put.
		## Numeric.put(b,ids,Numeric.take(b,ids)+bb)
                tools.vector.putAdd(b, ids, bb)

		
        ## explicit
        if self.weight.has_key('explicit'):
            
            oldArrayId1 = Numeric.take(oldArray, id1)
            oldArrayId2 = Numeric.take(oldArray, id2)
            
            tools.vector.putAdd(b, id1, -(self.explicit['cell 1 diag'][:self.interiorN] * oldArrayId1[:] + self.explicit['cell 1 offdiag'][:self.interiorN] * oldArrayId2[:]))
            tools.vector.putAdd(b, id2, -(self.explicit['cell 2 diag'][:self.interiorN] * oldArrayId2[:] + self.explicit['cell 2 offdiag'][:self.interiorN] * oldArrayId1[:]))

##            for i in range(self.interiorN):

##		b[id1[i]] -= self.explicit['cell 1 diag'][i] * oldArray[id1[i]] + self.explicit['cell 1 offdiag'][i] * oldArray[id2[i]]
##		b[id2[i]] -= self.explicit['cell 2 diag'][i] * oldArray[id2[i]] + self.explicit['cell 2 offdiag'][i] * oldArray[id1[i]]

            for boundaryCondition in self.boundaryConditions:

                LL,bb,ids = boundaryCondition.getContribution(self.explicit['cell 1 diag'],self.explicit['cell 1 offdiag'])
                oldArrayIds = Numeric.take(oldArray, ids)
                tools.vector.putAdd(b, ids, -LL * oldArrayIds)
                tools.vector.putAdd(b, ids, bb)
