#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "independentSourceTerm.py"
 #                                    created: 11/28/03 {11:36:25 AM} 
 #                                last update: 12/7/04 {2:57:30 PM} 
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
 #  2003-11-12 JEG 1.0 original
 # ###################################################################
 ##


from fipy.terms.sourceTerm import SourceTerm

class IndependentSourceTerm(SourceTerm):
    """
    Source term that does not depend on the solution variable. 
    Subtracted from the b vector.
    """
    def __init__(self, sourceCoeff, mesh):
        weight = {
	    'b vector': -1, 
	    'new value': 0, 
	    'old value': 0, 
	    'diagonal' : 0
	}
	
	SourceTerm.__init__(self, sourceCoeff, weight, mesh) 
