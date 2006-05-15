#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "tools.py"
 #                                    created: 11/17/03 {5:05:47 PM} 
 #                                last update: 5/15/06 {3:57:58 PM} 
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
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
 #  2003-11-17 JEG 1.0 original
 # ###################################################################
 ##

"""Vector utility functions that are inexplicably absent from Numeric
"""

import Numeric

from fipy.tools.dimensions.physicalField import PhysicalField
from fipy.tools import numerix

def crossProd(v1,v2):
    """Return vector cross-product of v1 and v2.
    """
    return PhysicalField(value = [v1[1] * v2[2] - v1[2] * v2[1],
			    v1[2] * v2[0] - v1[0] * v2[2],
			    v1[0] * v2[1] - v1[1] * v2[0]])
	
def sqrtDot(v1,v2):
    """Return square root of vector dot-product of v1 and v2.
    
	Usually used with v1==v2 to return magnitude of v1.
    """
#     return Numeric.sqrt(Numeric.dot(v2,v1))
##     return Numeric.sqrt(v1.dot(v2))
    ## We can't use Numeric.dot on quantities with units
##     return Numeric.sqrt(Numeric.sum(v1*v2))
    return numerix.sqrt(abs(numerix.sum(v1 * v2)))

def _putAddPy(vector, ids, additionVector, mask = None):
    additionVector = Numeric.array(additionVector)
    if mask is None:
        for i in range(len(ids)):
            vector[ids[i]] += additionVector[i]
    else:
        for i in range(len(ids)):
            if not mask[i]:
                vector[ids[i]] += additionVector[i]

def _putAddIn(vector, ids, additionVector):
    from fipy.tools.inline import inline
    inline._runInlineLoop1("""
	vector(ids(i)) += additionVector(i);
    """,
    vector = vector, ids = ids, additionVector = Numeric.array(additionVector),
    ni = len(ids))

def putAdd(vector, ids, additionVector):
    """ This is a temporary replacement for Numeric.put as it was not doing
    what we thought it was doing.
    """
    from fipy.tools.inline import inline
    inline._optionalInline(_putAddIn, _putAddPy, vector, ids, additionVector)

def prune(array, shift, start = 0):
    """
    removes elements with indices i = start + shift * n
    where n = 0, 1, 2, ...

        >>> prune(numerix.arange(10), 3, 5)
        [0,1,2,3,4,5,6,7,8,9,]
        >>> prune(numerix.arange(10), 3, 2)
        [0,1,3,4,6,7,9,]
        >>> prune(numerix.arange(10), 3)
        [1,2,4,5,7,8,]
        >>> prune(numerix.arange(4, 7), 3)
        [5,6,]

    """

    takeArray = numerix.nonzero(numerix.arange(len(array)) % shift != start)
    return Numeric.take(array, takeArray)

def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__":
    _test() 
