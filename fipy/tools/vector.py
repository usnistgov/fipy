#!/usr/bin/env python

## 
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "tools.py"
 #                                    created: 11/17/03 {5:05:47 PM} 
 #                                last update: 12/19/03 {4:08:07 PM} 
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

def crossProd(v1,v2):
    """Return vector cross-product of v1 and v2.
    """
    return Numeric.array([v1[1] * v2[2] - v1[2] * v2[1],
			    v1[2] * v2[0] - v1[0] * v2[2],
			    v1[0] * v2[1] - v1[1] * v2[0]],'d')
	
def sqrtDot(v1,v2):
    """Return square root of vector dot-product of v1 and v2.
    
	Usually used with v1==v2 to return magnitude of v1.
    """
    return Numeric.sqrt(Numeric.dot(v2,v1))

def arraySqrtDot(a1,a2):
    """Return array of square roots of vector dot-products
    for arrays a1 and a2 of vectors v1 and v2
    
    Usually used with v1==v2 to return magnitude of v1.
    """
    ## We can't use Numeric.dot on an array of vectors
    return Numeric.sqrt(Numeric.sum((a1*a2)[:],1))

def putAdd(vector, ids, additionVector):
    """ This is a temporary replacement for Numeric.put as it was not doing
    what we thought it was doing.
    """

    for i in range(len(ids)):
        vector[ids[i]] += additionVector[i]
    
    
