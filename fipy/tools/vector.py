#!/usr/bin/env python

## 
 # ###################################################################
 #  PyFiVol - Python-based finite volume PDE solver
 # 
 #  FILE: "tools.py"
 #                                    created: 11/17/03 {5:05:47 PM} 
 #                                last update: 2/2/04 {11:31:44 AM} 
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

"""Vector utility functions that are inexplicably absent from Numeric
"""

import Numeric

from fivol.tools.dimensions.physicalField import PhysicalField
import fivol.tools.array

import fivol.inline.inline as inline

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
    return fivol.tools.array.sqrt(fivol.tools.array.sum(v1*v2))

def arraySqrtDot(a1, a2, result):
    """Return array of square roots of vector dot-products
    for arrays a1 and a2 of vectors v1 and v2
    
    Usually used with v1==v2 to return magnitude of v1.
    """
    ## We can't use Numeric.dot on an array of vectors
##     return Numeric.sqrt(Numeric.sum((a1*a2)[:],1))
##    return fivol.tools.array.sqrt(fivol.tools.array.sum((a1*a2)[:],1))
    return inline.optionalInline(_arraySqrtDotIn, _arraySqrtDotPy, a1, a2, result)

def _arraySqrtDotPy(a1, a2, result):
    return fivol.tools.array.sqrt(fivol.tools.array.sum((a1*a2)[:],1))

def _arraySqrtDotIn(a1, a2, result):
    ni, nj = Numeric.shape(a1)
    inline.runInlineLoop1("""
        int j;
        result(i) = 0.;
        for (j = 0; j < nj; j++)
        {
            result(i) += a1(i,j) * a2(i,j);
        }
        result(i) = sqrt(result(i));
    """,result = result, a1 = a1, a2 = a2, ni = ni, nj = nj) 
    return result

def _putAddPy(vector, ids, additionVector):
    for i in range(len(ids)):
	vector[ids[i]] += additionVector[i]

def _putAddIn(vector, ids, additionVector):
    inline.runInlineLoop1("""
	vector(ids(i)) += additionVector(i);
    """, 
    vector = vector, ids = ids, additionVector = additionVector,
    ni = len(ids))

def putAdd(vector, ids, additionVector):
    """ This is a temporary replacement for Numeric.put as it was not doing
    what we thought it was doing.
    """
    inline.optionalInline(_putAddIn, _putAddPy, vector, ids, additionVector)

def prune(array, shift, start = 0):
    """
    removes elements with indices i = start + shift * n
    where n = 0, 1, 2, ...
    """
    
    takeArray = [x for x in range(len(array)) if (x % shift) != start]
    return Numeric.take(array, takeArray)
