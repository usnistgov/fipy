#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "explicitDiffusionTerm.py"
 #                                    created: 11/27/03 {11:39:03 AM} 
 #                                last update: 12/7/04 {2:48:20 PM} 
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
 #  2003-11-13 JEG 1.0 original
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

from fipy.terms.diffusionTerm import DiffusionTerm

class ExplicitDiffusionTerm(DiffusionTerm):
    r"""
    The discretization for the `ExplicitDiffusionTerm` is given by

    .. raw:: latex

       $$ \int_V \nabla \cdot (\Gamma\nabla\phi) dV \simeq \sum_f \Gamma_f
       \frac{\phi_A^\text{old}-\phi_P^\text{old}}{d_{AP}} A_f $$ where $\phi_A^\text{old}$ and
       $\phi_P^\text{old}$ are the old values of the variable. The term is
       added to the RHS vector and makes no contribution to
       the solution matrix.
    """
    
    def _getWeight(self, mesh):
	return {
	    'explicit':{
		'cell 1 diag':    -1, 
		'cell 1 offdiag':  1, 
		'cell 2 diag':    -1, 
		'cell 2 offdiag':  1
	    }
	}
	
	


